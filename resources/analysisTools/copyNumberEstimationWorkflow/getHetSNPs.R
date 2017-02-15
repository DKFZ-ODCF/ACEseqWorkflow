#!/usr/bin
#kortine Kleinheinz k.kleinheinz@dkfz.de
#first draft to estimate control allele frequencies based on the tumor BAF

library(getopt)
library(data.table)

script_dir = dirname(get_Rscript_filename())
source(paste0(script_dir,"/qq.R"))
source(paste0(script_dir, "/getopt.R"))


readfile <- function(filename, chrom){
		colNamesAllele <- c("CHROM", "POS", "An", "Bn", "At", "Bt", "haplo")
		write( paste("Reading ", chrom," from ", file_snp, "...\n" ), stderr())
		alleleList.chr <- try( fread( paste( "tabix ", file_snp, chrom ), header=FALSE, sep='\t', verbose=FALSE, showProgress=FALSE )  )
		if ( is.data.frame(alleleList.chr)  ){
			colnames( alleleList.chr ) = colNamesAllele
			alleleList.chr$CHROM <- gsub("^chr", "", alleleList.chr$CHROM)
			return(alleleList.chr)
		}else{
			write( paste(chrom," not found in ", file_snp, "skipping!\n"), stderr())
			return(NULL)
		}
}

writetable <- function(data, newFile, write.head){

  colNamesAllele <- c("#chr", "startPos", "Anormal", "Bnormal", "Atumor", "Btumor", "haplotype", "genotype")
  colnames(data) <- colNamesAllele
  format = format(data, scientific = FALSE, trim = TRUE)

  write.table(format, file = "", sep='\t', col.names=write.head,  row.names=FALSE, quote=FALSE)
}

estimateControlBAF <- function(data, minHet, maxHet){
	# remove uncovered regions
	selRm <- which(data$At == 0 & data$Bt ==0)
	if (length(selRm) >0)
		data <- data[-selRm,]

	tBAF <- data$Bt/(data$At + data$Bt)
	data$An <- 15
	data$Bn <- 15
	data$genotype <- "0/1"

	selHomoA <- which(tBAF < minHet)
	selHomoB <- which(tBAF > maxHet)

	data$An[selHomoA] <- 30
	data$Bn[selHomoA] <- 0
	data$genotype[selHomoA] <- "0/0"

	data$An[selHomoB] <- 0
	data$Bn[selHomoB] <- 30
	data$genotype[selHomoB] <- "1/1"

	return(data)
}

hetMin <- 0.1
hetMax <- 0.9

getopt2( matrix(c( 'file_snp',           's', 1, "character", "input coverage data, dbSNP",               # input, /ibios/co02/bludau/ACEseq/medullo_pediatric/MBBL8/all.snp.tab
		   'hetMin',	         'l','2',"numeric", "minimum BAF in tumor to be heterozygous",
		   'hetMax',	         'h','2',"numeric", "maximum BAF in tumor to be hetereozygous"
                ), ncol = 5, byrow = TRUE))
      
write(qq("file_snp: @{file_snp}\n\n"), stderr())
write("\n", stderr())



chromosomes= 1:24
write.head <<- TRUE
#print header
#cat(paste( c("#chr", "startPos", "Anormal", "Bnormal", "Atumor", "Btumor", "haplotype", "genotype"), sep="\t" ), "\n" )

tmp <-sapply( chromosomes, function(chr){
	allele <-  readfile(file_snp, chr)
	if( ! is.null(allele)){
		allele <-  estimateControlBAF(allele, hetMin, hetMax)
		writetable( allele, newFile = "|" , write.head=write.head)
		write.head <<-FALSE
	}
	
})
