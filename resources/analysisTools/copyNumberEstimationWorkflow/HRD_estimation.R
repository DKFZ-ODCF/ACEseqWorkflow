#!/usr/bin/R

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the MIT License (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

################################################################################
## Annotate segments with DUP/DEL/LOH/TCNneutral/NA based on full ploidy and CN 
## segment 
################################################################################

args <- commandArgs(trailingOnly = TRUE)
segmentfile <- args[1]
mergedfile <- args[2]
patientsex <- args[3]
ploidy <- as.integer(args[4])
tcc <- as.numeric(args[5])
pid <- args[6]
outfile <- args[7]
contributingSegmentsFile <- args[8]
centromerFile <- args[9]
cytobandsFile <- args[10]
pipelineDir <- args[11]
if( length(args)>11 ){
	cutoff <- as.numeric(args[12])
}else{
	cutoff <- 0.7
}


# define subtelomericCytobands which contains telomeric regions as they are defined in TelomereHunter.
cytoband.df <- read.table(cytobandsFile, header=F, stringsAsFactor=FALSE)
colnames(cytoband.df) = c("chrom","start", "end","cytoband","giemsa")
cytoband.df$chrom = gsub("chr", "", cytoband.df$chrom)

#cytoband.df$chrom = gsub("X", 23, cytoband.df$chrom)
#cytoband.df$chrom = gsub("Y", 24, cytoband.df$chrom)
cytoband.df$chrom = as.character(cytoband.df$chrom)
cytoband.df = cytoband.df[order(cytoband.df$chrom, cytoband.df$start),]

cytobandPerChr <- split(cytoband.df, cytoband.df$chrom)
subtelomericCytobands = lapply(cytobandPerChr, function(currentBands) {
	return(rbind(currentBands[1,],currentBands[nrow(currentBands),]))
})



source( file.path(pipelineDir, "annotateCNA.R") )

segments.df <- read.table(segmentfile, header=TRUE, stringsAsFactor=FALSE)
centromers <- read.table(centromerFile, header=FALSE, sep="\t", stringsAsFactor=FALSE)
colnames(centromers) <- c("chromosome", "start","end", "arm", "cytoband")



#calculate aberrant fractions
totallength <- sum( as.numeric(segments.df$length) )
totalAberrant <- sum( as.numeric(segments.df$length[segments.df$CNA.type != "TCNneutral"]),
			 na.rm=TRUE)
totalLost <- sum( as.numeric(segments.df$length[grep( "(DEL)|(HomeDel)",segments.df$CNA.type)]) )
totalLostLOH <- sum( as.numeric(segments.df$length[ which( grepl( "(DEL)|(HomeDel)|(LOH)",segments.df$CNA.type) & ! grepl("DUP", segments.df$CNA.type ) ) ] ) )
totalLOH <- sum( as.numeric(segments.df$length[ which( grepl( "LOH",segments.df$CNA.type) & ! grepl("DUP", segments.df$CNA.type ) ) ] ) )
totalGain <- sum( as.numeric(segments.df$length[grep( "DUP",segments.df$CNA.type)]) )

fractionAberrant <- totalAberrant/totallength
fractionGain <- totalGain/totallength
fractionLoss <- totalLost/totallength
fractionLossLOH <- totalLostLOH/totallength
fractionLOH <- totalLOH/totallength

merged.df <- read.table( mergedfile, header=TRUE, sep="\t" )
merged.df <- annotateCNA( seg.df = merged.df, ploidy=ploidy, cut.off = cutoff,
                            TCN.colname = "tcnMean", c1Mean.colname = "c1Mean",
                            c2Mean.colname = "c2Mean", sex=patientsex )

#- HRD score: 15MB
merged.df$roundTCN <- round(merged.df$tcnMean)

tcnStatePerChrom <- sapply(unique(merged.df$chromosome), function(i){
		sel <- which(merged.df$chromosome==i)
		nrow(unique( merged.df[sel, c("CNA.type", "roundTCN")]))
})
names(tcnStatePerChrom) <- unique(merged.df$chromosome)
selNoChangeChr <- names(tcnStatePerChrom)[which(tcnStatePerChrom==1)]

numberHRD        <- length( which(grepl("LOH", segments.df$CNA.type) & segments.df$length>15000000) )
numberHomoDel    <- length( which(grepl("HomoDel", segments.df$CNA.type) ) )

if(length(selNoChangeChr) != length(unique(merged.df$chromosome)) ){
  merged.df <- merged.df[! merged.df$chromosome %in% selNoChangeChr,]
  index.HRDSmooth = which(grepl("LOH", merged.df$CNA.type) & merged.df$length>15000000 )
  numberHRDSmooth  <- length( index.HRDSmooth )
  write.table( merged.df[index.HRDSmooth,], contributingSegmentsFile, sep="\t", row.names=FALSE, quote=FALSE )

  numberHRDLoss    <- length( which(grepl("(LOH)|(DEL)", merged.df$CNA.type) & 
                                   merged.df$length>15000000 ) )

  TAI=0
  segmentsPerChr <- split(merged.df, merged.df$chromosome)
  TAI <- sum( sapply( segmentsPerChr, function(segs){
		TAI_chr = 0
		segs$A = as.character(segs$A)
		segs$B = as.character(segs$B)
		# minimum length of TAI candidate regions: 11 Mbp
		# (according to Suppl. Doc1; Timms et al., Myriad Genetics Inc., Association of BRCA1/2 Defects with Genomic Scores Predictive of DNA Damage Repair Deficiency Among Breast Cancer Subtypes
		currentChrom = segs$chromosome[1]
  		endOfpArmSubtelomericRegion = subtelomericCytobands[[currentChrom]][1,"end"]
  		startOfqArmSubtelomericRegion = subtelomericCytobands[[currentChrom]][2,"start"]
		if( # segment starts in telomeric region
			segs[1,"start"] < endOfpArmSubtelomericRegion &
			# segment does not cross centromere
			segs[1,"end"] < centromers$end[centromers$chromosome == currentChrom][1] &
			segs[1,"length"] >= 11000000 &
			# imbalanced ?
			!is.na(segs[1,"A"]) & !is.na(segs[1,"B"]) &
			segs[1,"A"] != segs[1,"B"] ) {
					TAI_chr <- TAI_chr + 1
		}
		if( # segment starts after centromeric region
			segs[nrow(segs),"start"] > centromers$start[centromers$chromosome == currentChrom][2] &
			# segment ends within subtelomeric region
		    segs[nrow(segs),"end"] > startOfqArmSubtelomericRegion &
			segs$length[nrow(segs)] >= 11000000 &
			# imbalanced ?
			!is.na(segs[nrow(segs),"A"]) & !is.na(segs[nrow(segs),"B"]) &
			segs[nrow(segs),"A"] != segs[nrow(segs),"B"] ) {
					TAI_chr <- TAI_chr + 1
		}
		return(TAI_chr)
  }) )

  # LST score
  i=1
  LST=0
  for( j in 2:nrow(merged.df) ){
  	if( merged.df$chromosome[i] != merged.df$chromosome[j] | merged.df$length[i] < 10e6 | merged.df$length[j] <10e6 ){
  		i=i+1
  		next
  	}
  	if( merged.df$start[j] - merged.df$end[i] > 1){
  		i=i+1
  		next
  	}
  	if( abs(merged.df$tcnMean[j] - merged.df$tcnMean[i]) > cutoff){
  		LST=LST+1
  	}
  	i=i+1
  }
} else {
	# all chromosomes have only 1 state
	numberHRDSmooth <- 0
	numberHRDLoss <- 0
	LST <- 0
	TAI=0
}

out.data <- data.frame( pid, fractionAberrant, fractionGain, fractionLoss, 
                        fractionLossLOH, fractionLOH, numberHRDSmooth, numberHRD,
                        numberHomoDel, LST, numberHRDLoss, TAI )
write.table( out.data, outfile, sep="\t", row.names=FALSE, quote=FALSE )
