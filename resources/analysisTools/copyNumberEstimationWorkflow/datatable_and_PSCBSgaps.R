#!/usr/bin/R

library(getopt)

script_dir = dirname(get_Rscript_filename())
#source(paste0(script_dir,"/qq.R"))
#source(paste0(script_dir, "/getopt.R"))

################################################################################
## I wrap getopt into a function, which can check options, default values
## and export options to global environment
################################################################################
wd = getwd()
# set default values 
file_beta          = paste0(wd,"/densityBeta.pdf")
file_knownSegments = paste0(wd,"/knownSegments.txt")
file_data          = paste0(wd,"/pscbs_data.txt")
min_gap_length     = 500000
libloc=""
# variables in the first column will be exported into global environment
# I add the fifth column because it can explain the options in the first column
spec <- matrix( c('file_cnv',           'c', 1, "character", #"input coverage data, 1 kb windows", 
                 'file_snp',           's', 1, "character", #"input coverage data, dbSNP",               # input
		 'file_sex',	       'g','1',"character", #"file with sex of patient",		# input 
                 'file_beta',          'p', 2, "character", #"control plots path, should be .pdf file",  # output
                 'file_knownSegments', 'k', 2, "character", #"output file for known segments",           # output
                 'file_data',          'd', 2, "character", #"output file for SNPS",                     # output
                 'min_gap_length',     'l', 2, "numeric",   #"minimal gap length, default is 500000",
		 'libloc',		'o',2, "character" #"location of PSCBS package"
                ), ncol = 4, byrow = TRUE)
      
opt = getopt(spec);
for(item in names(opt)){
       assign( item, opt[[item]])
}

cat(paste0("file_cnv:  ",file_cnv, "\n\n"))
cat(paste0("file_snp:  ",file_snp, "\n\n"))
cat(paste0("file_beta: ",file_beta, "\n\n"))
cat(paste0("file_knownSegments: ", file_knownSegments, "\n\n"))
cat(paste0("file_data: ", file_data, "\n\n"))
cat(paste0("file_sex: ", file_sex, "\n\n"))
cat(paste0("min_gap_length: ", min_gap_length, "\n\n"))
cat("\n")

if ( libloc == "" | libloc == TRUE ){
	libloc = NULL
}

library(PSCBS, lib.loc=libloc)

#input coverage data, 1 kb windows
cat(paste0("reading ", file_cnv, "...\n\n"))
sample = read.table(file_cnv, sep = "\t",as.is = TRUE, header=T)
          
#input coverage data, dbSNP
cat(paste0("reading ",file_snp, "...\n\n"))
colNamesAllele = c("chr", "pos", "Anormal1", "Bnormal1", "Atumor1", "Btumor1", "haplotype")

#if patient is female remove all Y chromosome windows (will lead to exclusion of SNPs during merge) 
sex = read.table(file_sex, header=FALSE, stringsAsFactors=FALSE)[,1]

chromosomes = 1:24

if (sex == "female") {
	sample <- sample[sample$chromosome != 'chr24' & sample$chromosome != 24,]
	chromosomes = 1:23
}

allele = lapply(chromosomes, function(chr){
		chrom= paste0(chr)
		cat( "Reading ", chrom," from ", file_snp, "...\n" )
		alleleList.chr <- try( read.table( pipe( paste( "tabix ", file_snp, chrom ) ), header=FALSE, sep='\t' )  )
		if ( is.data.frame(alleleList.chr)  ){
			colnames( alleleList.chr ) = colNamesAllele
			alleleList.chr
		}else{
			cat(chr," not found in ", file_snp, "skipping!\n")
			NULL
		}
	})

allele <- do.call(rbind, allele)



#calculate total number of reads for normalization
sumSampleN = sum(as.numeric(sample$normal))
sumSampleT = sum(as.numeric(sample$tumor))

#B-allele frequency (NOT normalized)   
allele$BAFnormal = allele$Bnormal1 / (allele$Anormal1 + allele$Bnormal1)       
allele$BAFtumor = allele$Btumor1 / (allele$Atumor1 + allele$Btumor1)
          
#SNP-position rounded to window start (1kb bins)
## 1kb bin, hard coded
#allele$covWindow = round(allele$pos, digits= -4) #bludau script: rounded to lower 10000
allele$covWindow = floor(allele$pos/10000)*10000+1          
sample$start = floor(sample$start/10000)*10000+1          
#SNP values
tableallele = data.frame(a = allele$covWindow,
                         betaT = allele$BAFtumor, 
                         betaN = allele$BAFnormal, 
                         x = allele$pos, 
                         chromosome = allele$chr, 
                         Atumor = allele$Atumor1, 
                         Btumor = allele$Btumor1, 
                         Anormal = allele$Anormal1, 
                         Bnormal = allele$Bnormal1,
			 haplotype = allele$haplotype)


#tablesample = data.frame(a = round(sample$startPos, digits = -4), #bludau: already rounded in previous script
tablesample = data.frame(a = sample$start,
                         chromosome = sample$chromosome,
                         CT = sample$covR,
                         covT = sample$tumor,
                         covN = sample$normal,
			 covR_raw = sample$covR_raw)
        
#remove "chr" in front of the chromosome number
tablesample$chromosome = gsub("^chr", "", tablesample$chromosome)
tableallele$chromosome = gsub("^chr", "", tableallele$chromosome)

cat(paste0( nrow(tableallele), " rows for tableallele.\n\n"))
cat(paste0( nrow(tablesample), " rows for tablesample\n\n"))

################################################################################
## here I used sqldf instead of data.table to run SQL command to merge two data 
## frames. Because data.table seems also using SQLite to do data process but it 
## is not strong enough as sqldf to utilize all functionality of SQLite
##
## I checked that output by sqldf and data.table are all same.
################################################################################
cat("merging data frames...\n")
completeTable <- merge(  x=tableallele,
			 y=tablesample,
			 by=c('a', 'chromosome' ),
			 sort=FALSE )

completeTable <- completeTable[order(completeTable$chromosome, completeTable$x),] 

cat("generating control plots.\n")       
 
#control plots of betaT and betaN 
################################################################################
## I put two density plot in one figure. Also the format of figure output is
## changed. I prefer PDF format but PNG can also be generated according to 
## the suffix name of figure
################################################################################


require(ggplot2)
p <- ggplot( data=completeTable, aes(x=betaT, y=..scaled..), colour='red') + geom_density()
p <- p + geom_density(aes(x=betaN, y=..scaled..), colour='blue')
p <- p + scale_colour_manual( name=c("Tumor","Normal") , breaks=c('betaT', 'betaN'), values=c('red','blue') )
p <- p + ggtitle("control plots of betaT and betaN") + ylab("scaled density")
ggplot2::ggsave(file_beta, p, width=6, height=6, units='cm')

#PSCBS
cat("PSCBS\n")

data = dropSegmentationOutliers(completeTable) #removes single outliers
gaps = findLargeGaps(data, minLength = min_gap_length) #regions larger than 5 Mb that do not contain enough data points (SNPs) are detected (e.g. centromeres)
knownSegments = gapsToSegments(gaps) #the gaps are stored as predefined segments in "knownSegments)

write.table(knownSegments, file = file_knownSegments, sep = "\t", row.names = FALSE, quote = FALSE )

write.table(data, file = pipe( paste0("bgzip >", file_data) ), sep = "\t", row.names = FALSE, quote = FALSE )
