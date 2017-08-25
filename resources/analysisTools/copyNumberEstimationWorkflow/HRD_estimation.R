#!/usr/bin/R

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
centromerFile <- args[8]
pipelineDir <- args[9]
if( length(args[10])>9 ){
	cutoff <- as.numeric(args[10])
}else{
	cutoff <- 0.7
}

source( file.path(pipelineDir, "annotateCNA.R") )

segments.df <- read.table(segmentfile, header=TRUE, stringsAsFactor=FALSE)
centromers <- read.table(centromerFile, header=FALSE, sep="\t", stringsAsFactor=FALSE)
colnames(centromers) <- c("chromosome", "start","end", "arm", "cytoband")

#calculate aberrant fractions
totallength <- sum( as.numeric(segments.df$length) )
totalAberrant <- sum( as.numeric(segments.df$length[segments.df$CNA.type != "TCNneutral"]) )
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
		length(unique( merged.df[sel, c("CNA.type", "roundTCN")]))
		})
names(tcnStatePerChrom) <- unique(merged.df$chromosome)
selNoChangeChr <- names(tcnStatePerChrom)[which(tcnStatePerChrom==1)]

numberHRD        <- length( which(grepl("LOH", segments.df$CNA.type) & segments.df$length>15000000) )
numberHomoDel    <- length( which(grepl("HomoDel", segments.df$CNA.type) ) )

if(length(selNoChangeChr) != length(unique(merged.df$chromosome)) ){
  merged.df <- merged.df[! merged.df$chromosome %in% selNoChangeChr,]

  numberHRDSmooth  <- length( which(grepl("LOH", merged.df$CNA.type) & 
                             merged.df$length>15000000 ) )
  numberHRDLoss    <- length( which(grepl("(LOH)|(DEL)", merged.df$CNA.type) & 
                                   merged.df$length>15000000 ) )

  TAI=0
  segmentsPerChr <- split(merged.df, merged.df$chromosome)
  TAI <- sum( sapply( segmentsPerChr, function(segs){
		TAI_chr = 0

		if( segs$end[1] <
			centromers$end[centromers$chromosome == segs$chromosome[1]][1] &
		    segs$start[1] <
			centromers$end[centromers$chromosome == segs$chromosome[1]][1] ) {
			TAI_chr <- TAI_chr + sum( grepl( "(DEL)|(DUP)|(LOH)", segs$CNA.type[1]) )
		}
		if( segs$start[nrow(segs)] >
			centromers$start[centromers$chromosome == segs$chromosome[1]][2] & 
		    segs$end[nrow(segs)] >
			centromers$start[centromers$chromosome == segs$chromosome[1]][2] ) {
			TAI_chr <- TAI_chr + sum( grepl( "(DEL)|(DUP)|(LOH)", segs$CNA.type[nrow(segs)]) )
		}
		TAI_chr
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
}else{
  numberHRDSmooth <- 0
  numberHRDLoss <- 0
  LST <- 0
}

out.data <- data.frame( pid, fractionAberrant, fractionGain, fractionLoss, 
                        fractionLossLOH, fractionLOH, numberHRDSmooth, numberHRD,
                        numberHomoDel, LST, numberHRDLoss, TAI )
write.table( out.data, outfile, sep="\t", row.names=FALSE, quote=FALSE )
