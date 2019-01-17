#!/usr/bin/R

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
library(getopt)

script_dir = dirname(get_Rscript_filename())

wd = getwd()
libloc=NULL
# set default values 
file_fit  = paste0(wd, "/fit.txt")
nperm     = 200
min.width = 2000
undo.SD   = 8
trim      = 1e-6
alphaTCN  = 1e-6
alphaDH   = 1e-6
h         = 0
sv        = "true"
nocontrol = FALSE 
################################################################################
## description field should be filled here
################################################################################
spec <- matrix(c('file_data',        'd', 1, "character", #"", # input
                 'file_breakpoints', 'b', 1, "character", #"", # input
		         'chrLengthFile',    'x', 1, "character", #"chromosome length file",
                 'file_fit',         'f', 2, "character", #"", # output
                 'nperm',            'p', 2, "numeric"  , #  "",
                 'minwidth',         'w', 2, "integer"  , #  "",
                 'undo.SD',          's', 2, "integer"  , #  "",
                 'trim',             't', 2, "numeric"  , #  "",
                 'alphaTCN',         'T', 2, "numeric"  , #  "",
                 'alphaDH',          'D', 2, "numeric"  , #  "",
                 'h',                'h', 2, "numeric"  , #  "",
                 'sv',               'c', 0, "character", #"whether svs are used, should be specified",
                 'libloc',           'l', 2, "character", #"location of package",
		         'nocontrol',	     'n', 2, "logical"    #"is the sample run without control"
                ), ncol = 4, byrow = TRUE)
 
opt = getopt(spec);
for(item in names(opt)){
       assign( item, opt[[item]])
}
    
    
cat(paste0("file_data: ", file_data, "\n\n"))
cat(paste0("file_breakpoints: ", file_breakpoints, "\n\n"))
cat(paste0("chrLengthFile: ", chrLengthFile, "\n\n"))
cat(paste0("file_fit: ", file_fit, "\n\n"))
cat(paste0("nperm: ", nperm, "\n\n"))
cat(paste0("minwidth: ", minwidth, "\n\n"))
cat(paste0("undo.SD: ", undo.SD, "\n\n"))
cat(paste0("trim: ", trim, "\n\n"))
cat(paste0("alphaTCN: ", alphaTCN, "\n\n"))
cat(paste0("alphaDH: ", alphaDH, "\n\n"))
cat(paste0("h: ", h, "\n\n"))
cat(paste0("sv: ", sv, "\n\n"))
cat(paste0("libloc: ", libloc, "\n\n"))
cat(paste0("nocontrol: ", nocontrol, "\n\n"))
cat("\n")

if( libloc == ""  | libloc == TRUE )
    libloc=NULL

library(DNAcopy)
originalSegmentFunction <- get("segment", mode = "function", envir = getNamespace("DNAcopy"))
segment.CODE = deparse(originalSegmentFunction)
lineToModify = grep("stop\\(\"minimum segment width should be between 2 and 5\"\\)", segment.CODE)
segment.CODE[lineToModify] = "#HACK... - do not stop here"
modifiedSegmentFunction = eval(parse(text = segment.CODE))
R.utils::reassignInPackage("segment", "DNAcopy", modifiedSegmentFunction, keepOld=F)

library(PSCBS, lib.loc=libloc)
# read datatable                        
cat(paste0("reading ",file_data, "...\n\n"))
colNamesData <- c( "a", "chromosome", "betaT", "betaN", "x", "Atumor", "Btumor", "Anormal", "Bnormal", "haplotype", "CT", "covT", "covN" )
chromosomes = 1:24

data = lapply( chromosomes, function(chr){
		cat( "Reading chr ", chr," from ", file_data, "...\n" )
		dataList.chr <- try( read.table( pipe( paste( "tabix", file_data, chr ) ), header=FALSE, sep='\t', as.is=TRUE ) )
		if ( is.data.frame(dataList.chr)  ){
			colnames( dataList.chr ) = colNamesData
			if (nocontrol)
				dataList.chr$betaN <- jitter(dataList.chr$betaN)
			dataList.chr
		}else{
			cat(chr," not found in ", file_data, "skipping!\n")
			NULL
		}
	})

data <- do.call(rbind, data)

chrLength = read.table(chrLengthFile, header = FALSE, as.is = TRUE)
chrLength = data.frame(chrLength)
chrLength$V1 <- gsub('chr','',chrLength$V1)

sel = which(chrLength$V1 == "X")
chrLength$V1[sel] = 23
sel = which(chrLength$V1 == "Y")
chrLength$V1[sel] = 24

# read knownSegments
################################################################################
## I think `sv` should be a logical variable.
################################################################################
cat(paste0("reading ", file_breakpoints, "...\n\n"))
knownSegments = read.table(file_breakpoints, sep = "\t", as.is = TRUE, header=FALSE ) 
colnames(knownSegments) = c("chromosome", "start", "end", "length")

    
# start pscbs
cat("start pscbs\n")
#library(PSCBS, lib.loc="/home/bludau/R/x86_64-unknown-linux-gnu-library/2.13")
# lib.loc needs to be specified since source code has been modified
# remove a line in `DNAcopy::segment`
#body(segment)[[4]] = substitute(print("hack")) # minimum width 


fit = segmentByPairedPSCBS(data, knownSegments = knownSegments, tbn = FALSE, 
        nperm = nperm, min.width = min.width, undoTCN = 5, undoDH = 5, 
        undo.splits = "sdundo", undo.SD = undo.SD, trim = trim, alphaTCN = alphaTCN, 
        alphaDH = alphaDH, verbose = FALSE, seed = 25041988)
    
if (h != 0) {
        fit = pruneByHClust(fit, h = h, merge = TRUE, update = TRUE, verbose = -10)
}
cat("Finished segmentByPairedPSCBS")

segments = getSegments(fit, simplify = TRUE)
cat("Finished getSegments")

#round coordinates with .5 to closest integer (down for end, up dor start)
if( any( segments$start == -Inf, na.rm=TRUE ) ){
	selStart <- which( segments$start == -Inf )
	segments$start[ selStart ] = 1
}

if( any( segments$end == Inf, na.rm=TRUE ) ){
	sel <- which( segments$end == Inf )
	segments$end[ sel ] <-  chrLength$V2[  sapply(sel, function(i) which( chrLength$V1  == segments$chromosome[i] ) )]
}

segments$start	<- as.integer(ceiling(segments$start))
segments$end	<- as.integer(floor(segments$end))

segments = format(segments, scientific = FALSE, trim = TRUE)
colnames(segments) <- gsub("^chromosome", "#chromosome", colnames(segments))

write.table(segments, file = file_fit, sep = "\t", row.names = FALSE, quote = FALSE )
