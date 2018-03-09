#!/usr/bin/R

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

library(getopt)
library(data.table)

################################################################################
## I wrap getopt into a function, which can check options, default values
## and export options to global environment
################################################################################

# variables in the first column will be exported into global environment
# I add the fifth column because it can explain the options in the first column
spec <- matrix(c('file_snp',           's', 1, "character" , #"input coverage data, dbSNP",              # input,
		 'file_sex',	       'g','1',"character", #"file with sex of patient",		# input 
		 'chrLengthFile',      'f','1',"character", #"file with sex of patient",		# input 
		 'pid',                'p','1',"character", #"patient identifier",			# input 
                 'plot_Dir',          'd', 2, "character"  #"output file for SNPS"                     # output
                ), ncol = 4, byrow = TRUE)

opt = getopt(spec);
for(item in names(opt)){
       assign( item, opt[[item]])
}
     
cat(paste0("file_snp: ",file_snp, "\n\n"))
cat(paste0("plot_Dir: ",plot_Dir, "\n\n"))
cat(paste0("file_sex: ",file_sex, "\n\n"))
cat("\n")

plotCoverage <- function(coverageTab, chromosomeBorders=NULL, chr=NULL, ylims=c(0,1.2) ) {
    
    labelposition=NULL
    if ( ! is.null(chromosomeBorders) ){
	  labelPosition <- (diff(chromosomeBorders)/2+chromosomeBorders[1:length(chromosomes)])/1e6
	  #p <- ggplot(data=coverageTab[seq(1, nrow(coverageTab),7),], aes(x=pos.adjusted/1e6, y=BAFnormal)) + geom_point(alpha=0.7) +
	   #   xlim(c(0,max(chromosomeBorders/1e6))) +ylim(ylims) +ylab('control BAF') + xlab('genomic coordinate (MB)')
	  #p <- p + geom_hline(yintercept=0.5)
	  plot(coverageTab$pos.adjusted/1e6, coverageTab$BAFnormal, xlab='genomic coordinate (MB)', 
            ylab="control BAF", main=paste0( "genome ", chr, "\n", "control BAF"), ylim=ylims, 
            xlim=c(0, max(chromosomeBorders)/1e6 ), pch=16, cex= 0.4, col=rgb(0,0,0,0.5) )
	    abline(v=chromosomeBorders/1e6)
	    text( labelPosition, replicate(length(labelPosition), 1.2 ), gsub("^chr", "", c(1:22,"X","Y")[1:length(chromosomes)] ) )
    } else{
	  plot(coverageTab$pos/1e6, coverageTab$BAFnormal, xlab='genomic coordinate (MB)', 
        	 ylab="control BAF", main=paste0( "chromosome", chr, "\n", "control BAF"), ylim=ylims, 
	         xlim=c(0, max(coverageTab$pos)/1e6 ), pch=16, cex= 0.4, col=rgb(0,0,0,0.5) )
#          if(chr==23)
#		abline(v=c(113570000/1e6,113750000/1e6))
		  
   }
    #add chromosome labels and borders
    #labels.data        <- data.frame(xPos=labelPosition )
    #labels.data$maxCov <- replicate(nrow(labels.data), 1.1)
    #labels.data$chr    <- chromosomes

    #vlines   <- geom_vline( xintercept= chromosomeBorders[2:length(chromosomeBorders)]/1e6, col="#000000", lwd = 0.2 )
    #p <- p + vlines + labs(x=NULL) + theme( axis.text.x=element_blank() )
    #p <- p + geom_text( data=labels.data,aes( x=xPos, y=maxCov, label=chr),size=10 ) 
    #p
        
}


adjustCoordinates <- function(sortedLengthTab, coverageList){

    addToStart <- 0
    chromosomeBorders <- 0
    
    for (chr in sortedLengthTab$chromosome ){
      
      chrList <- coverageList[[chr]]
      if( nrow(chrList) > 0 )
	      chrList$pos.adjusted <- chrList$pos + addToStart
	
      addToStart <- addToStart + sortedLengthTab$length[sortedLengthTab$chromosome==chr] 
      chromosomeBorders <- c(chromosomeBorders, addToStart)
      coverageList[[chr]] <- chrList
    }
    
  return( list(coverageList=coverageList, chromosomeBorders=chromosomeBorders) )
}

#read and sort chromosome length file
chrLengthTab = read.table(chrLengthFile, header = FALSE, as.is = TRUE)
chrLengthTab = data.frame(chrLengthTab)
colnames(chrLengthTab)  <- c("chromosome", "length", "info")[1:dim(chrLengthTab)[2]]
chrLengthTab$chromosome <- gsub('chr','',chrLengthTab$chromosome)
chrLengthTab$chromosome <- gsub('X', 23, chrLengthTab$chromosome)
chrLengthTab$chromosome <- gsub('Y', 24, chrLengthTab$chromosome)
chrLengthTab$chromosome <- as.numeric(chrLengthTab$chromosome)

#input coverage data, dbSNP
cat(paste0("reading ",file_snp, "...\n\n"))
colNamesAllele = c("chromosome", "pos", "Anormal1", "Bnormal1", "Atumor1", "Btumor1", "haplotype")

#if patient is female remove all Y chromosome windows (will lead to exclusion of SNPs during merge) 
sex = read.table(file_sex, header=FALSE, stringsAsFactors=FALSE)[,1]

chromosomes = c(1:24)

if (sex == "female") {
	chromosomes = c(1:23)
}
#chromosomes <- gsub("^", "chr",chromosomes)
chrLengthTab		<- chrLengthTab[sapply(chromosomes, function(i) which(chrLengthTab$chromosome==i)),]

allele = lapply(chromosomes, function(chr){
		chrom= paste0(chr)
		cat( "Reading ", chrom," from ", file_snp, "...\n" )
		alleleList.chr <- try( fread( paste( "tabix ", file_snp, chrom ), header=FALSE, sep='\t' )  )
		if ( is.data.frame(alleleList.chr)  ){
			colnames( alleleList.chr ) = colNamesAllele
			alleleList.chr$chromosome <- gsub("^chr", "", alleleList.chr$chromosome)
			alleleList.chr <- alleleList.chr[seq(1, nrow(alleleList.chr), 10),]
			alleleList.chr$BAFnormal = alleleList.chr$Bnormal1 / (alleleList.chr$Anormal1 + alleleList.chr$Bnormal1)       
                        png(paste0(plot_Dir,"/control_", pid, "_", chr,"_BAF.png"), width=2000, height=1000, type="cairo")
			    plotCoverage( alleleList.chr, chr=chr )
			dev.off()

			alleleList.chr
		}else{
			cat(chr," not found in ", file_snp, "skipping!\n")
			NULL
		}
	})

names(allele) <- chromosomes

allele <- adjustCoordinates( chrLengthTab, allele )
chromosomeBorders <- allele$chromosomeBorders
gc()

allele <- allele$coverageList
allele <- rbindlist(allele)
gc()


png(paste0(plot_Dir,"/control_", pid, "_wholeGenome_BAF.png"), width=2000, height=1000, type='cairo')
     plotCoverage( allele, chromosomeBorders=chromosomeBorders )
dev.off()

