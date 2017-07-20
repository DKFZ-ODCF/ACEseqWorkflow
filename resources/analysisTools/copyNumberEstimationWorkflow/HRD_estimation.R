#!/usr/bin/R
suppressMessages(library(getopt))

script_dir = dirname(get_Rscript_filename())
source(file.path( script_dir, "qq.R") )
source(file.path( script_dir, "getopt.R") )

cutoff <- 0.7
################################################################################
## description field should be filled here
################################################################################
getopt2(matrix(c('segmentfile',      'f', 1, "character", "comb_pro_extra_file",
		 'mergedfile',       'm', 1, "character", "smoothed comb pro extra file",
                 'patientsex', 	     's', 1, "character", "patient sex (male/female)",
		 'ploidy',	     'p', 1, "integer", "full ploidy",
		 'tcc',	             't', 1, "numeric", "tumor cell content",
                 'pid',              'i', 1, "character", "patient identifier", 
                 'outfile',          'o', 1, "character", "outfile for parameters", 
                 'pipelineDir',      'u', 1, "character", "path to pscbs_plot_functions.R to load annotateCNV function", 
		 'cutoff',	     'c', 2, "numeric", "required deviation from full ploidy to be counted as aberrant"
                ), ncol = 5, byrow = TRUE));
print(outfile)
source( file.path(pipelineDir, "pscbs_plots_functions.R") )

segments.df <- read.table(segmentfile, header=TRUE)

#calculate aberrant fractions
totallength <- sum( as.numeric(segments.df$length) )
totalAberrant <- sum( as.numeric(segments.df$length[segments.df$CNA.type != "neutral"]) )
totalLost <- sum( as.numeric(segments.df$length[grep( "(DEL)|(HomeDel)",segments.df$type)]) )
totalLostLOH <- sum( as.numeric(segments.df$length[ which( grepl( "(DEL)|(HomeDel)|(LOH)",segments.df$CNA.type) & ! grepl("DUP", segments.df$CNA.type ) ) ] ) )
totalLOH <- sum( as.numeric(segments.df$length[ which( grepl( "LOH",segments.df$CNA.type) & ! grepl("DUP", segments.df$CNA.type ) ) ] ) )
totalGain <- sum( as.numeric(segments.df$length[grep( "DUP",segments.df$CNA.type)]) )

fractionAberrant <- totalAberrant/totallength
fractionGain <- totalGain/totallength
fractionLoss <- totalLost/totallength
fractionLossLOH <- totalLostLOH/totallength
fractionLOH <- totalLOH/totallength

#TODO
#- also merge over gaps

merged.df <- read.table( mergedfile, header=TRUE, sep="\t" )
merged.df <- annotateCNA( seg.df = merged.df, ploidy=ploidy, cut.off = 0.7,
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

numberHRD        <- length( which(grepl("LOH", segments.df$type) & segments.df$length>15000000) )
numberHomoDel    <- length( which(grepl("HomoDel", segments.df$type) ) )

if(length(selNoChangeChr) != length(unique(merged.df$chromosome)) ){
  merged.df <- merged.df[! merged.df$chromosome %in% selNoChangeChr,]

  numberHRDSmooth  <- length( which(grepl("LOH", merged.df$type) & 
                             merged.df$length>15000000 ) )
  numberHRDLoss    <- length( which(grepl("(LOH)|(DEL)", merged.df$type) & 
                                   merged.df$length>15000000 ) )



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
                        numberHomoDel, LST, numberHRDLoss )
write.table( out.data, outfile, sep="\t", row.names=FALSE )
