#!/usr/bin/Rscript

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

swapAlleles <- function(segments, data, chr, blockPre, blockPost){

  blockFile <- paste0( blockPre, chr, ".", blockPost)
  blocks <- read.table( blockFile, header=F)
  colnames(blocks) <- c('chr', 'start', 'end', 'length')
  
  starts  <- segments$start
  ends    <- segments$end
 
  #get heterozygous SNPs from normal tissue
  selSub  <- which(data$betaN>0.3 & data$betaN<0.7 & ! is.na(data$betaT), )
  dataSub <- data[ selSub, ]
    
  # loop over all segments on chromosome get median BAF dor all haplogroups, if median different from overall median, 
  # swap allele frequencies for all alleles within haploblock
  for ( j in seq_along(starts) ){
    
    # get median BAF for tumor within each haplotype block in segment
    selB <- which( ( blocks$end <= ends[j] & blocks$end >= starts[j] ) | ( blocks$start <= ends[j] & blocks$start >= starts[j] ) ) #blocks in segments
    if (length(selB) > 0 ){
      blockSub <- blocks[ selB, ]
      selDat <- which( dataSub$SNP >= min(blockSub$start) & dataSub$SNP <= max(blockSub$end) )
      
      swap <- c()
      for ( r in seq_len( nrow(blockSub) ) ) {

	#find SNPs that are lying within block or within segment boundary and block end
	s <- blockSub$start[r] 
	e <- blockSub$end[r] 
	if ( blockSub$start[r] < starts[j] ){
		s <- starts[j] 
	}
	if ( blockSub$end[r] > ends[j] ){
		e <- ends[j]
	}

	temp <- which( dataSub$SNP[selDat] >= s & dataSub$SNP[selDat] <= e ) #snps in block

        if ( length(temp) > 0 ) {
	         dataSub$adjusted[selDat[temp]] <- 1
			
                 BAFs <- median( dataSub$betaT[selDat[temp]] )

    	         #check median BAF for each haplotype Block and adjust Allele frequencies if not >0.5
      	         if (BAFs < 0.5){
                        swap <- c(swap, temp)
      	         }
        }        	 
      }
      if ( length(swap) > 0){
        # swap A and B allele for those blocks which differ from median (should be reversed?)
        tmpA  <- list(dataSub$Atumor[selDat[swap]], dataSub$Anormal[selDat[swap]])
        dataSub$Atumor[selDat[swap]]  <- dataSub$Btumor[selDat[swap]]
        dataSub$Anormal[selDat[swap]] <- dataSub$Bnormal[selDat[swap]]
        dataSub$Btumor[selDat[swap]]  <- tmpA[[1]]
        dataSub$Bnormal[selDat[swap]] <- tmpA[[2]]
     }
    }
  } #finish segment loop


  dataSub$betaT <- dataSub$Btumor/(dataSub$Btumor+dataSub$Atumor)
  dataSub$betaN <- dataSub$Bnormal/(dataSub$Bnormal+dataSub$Anormal)

  #replace adjusted lines in data.frame
  data[ selSub, ] <- dataSub
  #remove those that haven't been adjusted due to laking genotype information
  return(data)
}
