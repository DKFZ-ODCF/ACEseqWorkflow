#!/usr/bin/R

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

### Kortine Kleinheinz
### 28/02/2014 k.kleinheinz@dkfz.de
### Determine the gender of the patient with simple ratio of per base coverage Y chromosome over All chromosomes
library(getopt)

script_dir = dirname(get_Rscript_filename())

wd = getwd()
# set default values 
min_Y_ratio=0.12
min_X_ratio=0.8
################################################################################
## Get coverage of Y chromosome in normal tissue to determine gender 
################################################################################
spec <- matrix(c('file_dataY',      'd', 1, "character", #"CNV file for Y chromsome",
                 'file_dataX',      'f', 1, "character", #"CNV file for X chromosome",
		 'file_size',	    's', 1, "character", #"chromosome length file",
		 'cnv_files',	    'c', 1, "character", #"cnv file for all chromosomes (use '*' instead of chr number)",	
		 'file_out',	    'o', 1, "character", #"outfile with determined sex of patient",
                 'min_Y_ratio',	    'y', 2, "numeric"  , #"minimum ratio to be exceeded for male patients",
		 'min_X_ratio',     'x', 2, "numeric"    #"minimum ratio to be exceeded for female patients"
                ), ncol = 4, byrow = TRUE)

opt = getopt(spec);
for(item in names(opt)){
       assign( item, opt[[item]])
}
    
cat(paste0("file_dataY: ",file_dataY, "\n\n"))
cat(paste0("file_dataX: ",file_dataX, "\n\n"))
cat(paste0("file_size: ",file_size, "\n\n"))
cat(paste0("file_out: ",file_out, "\n\n"))
cat(paste0("cnv_files: ",cnv_files, "\n\n"))
cat(paste0("min_Y_ratio: ",min_Y_ratio, "\n\n" ))
cat(paste0("min_X_ratio: ",min_X_ratio, "\n\n" ))
cat("\n")

#read chromosome length file and cnv file for Y
dataX <- read.table( file_dataX, sep='\t', header=FALSE )
dataY <- read.table( file_dataY, sep='\t', header=FALSE )
size <- read.table( file_size, header=FALSE, as.is=TRUE )
size$V1 <- gsub('chr', '', size$V1)

colnames(dataX) <- c( 'chr', 'pos', 'end', 'normal', 'tumor' )
colnames(dataY) <- c( 'chr', 'pos', 'end', 'normal', 'tumor' )
colnames(size) <- c( 'chr', 'length' )

# sum of per base coverage over whole genome
index <- which(colnames(dataY)=="normal")
cmd <- paste0("awk '{s+=$", index,"} END{print s}' ", cnv_files)

covY <- sum(as.numeric(dataY$normal))
covX <- sum(as.numeric(dataX$normal))
covAll <- as.numeric(system(cmd, intern=TRUE))

lengthY   <- size[size$chr=='Y', 'length']
lengthX   <- size[size$chr=='X', 'length']
lengthAll <- sum( as.numeric(size$length) )

covMale <- (covY/lengthY)/(covAll/lengthAll)
covFemale <- (covX/lengthX)/(covAll/lengthAll)
cat( paste0("covMale: ",covMale, "\n covFemale:", covFemale, "\n\n" ) )

if ( covFemale >= min_X_ratio ){
  if (covMale >= min_Y_ratio ){
	sex='klinefelter'
  }else{
    	sex='female'
  }
}else{
	sex='male'
}

write.table(data.frame(sex), file_out, row.names=FALSE, col.names=FALSE, quote=FALSE)
