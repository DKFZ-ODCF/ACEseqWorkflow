#!/usr/bin/R

library(ks)
library(getopt)

script_dir = dirname(get_Rscript_filename())
source(paste0(script_dir,"/qq.R"))
source(paste0(script_dir, "/getopt.R"))
source(qq("@{script_dir}/functions.R"))

minLim=0.47
maxLim=0.53
minCoverage=20
getopt2(matrix(c('file',     	'f', 1, "character", "", # input  /ibios/co02/bludau/ACEseq/medullo_pediatric/MBBL8/all_seg_2.txt
		 'gender',	'g', 1, "character", "",
                 'segments',	's', 1, "character", "", # input  /ibios/co02/bludau/ACEseq/medullo_pediatric/MBBL8/clustered_and_pruned_and_normal.txt
                 'segOut',    	'o', 1, "character", "",  # output segments file
		 'out',		'u', 1, "character", "", #output dir
		 'minLim',	'i', 2, "numeric"  , 'minimum peak value to be within balanced range',
		 'maxLim',	'a', 2, "numeric" , "maximum peak value to be within balanced range",
		 'minCoverage',	'm', 2, "numeric" , "minimum coverage in control for a SNP to be considered"
                ), ncol = 5, byrow = TRUE));

cat(qq("file: @{file}\n\n"))
cat(qq("gender_file: @{gender}\n\n"))
cat(qq("segments: @{segments}\n\n"))
cat(qq("segOut: @{segOut}\n\n"))
cat(qq("out: @{out}\n\n"))
cat(qq("minLim: @{minLim}\n\n"))
cat(qq("maxLim: @{maxLim}\n\n"))
cat(qq("minCoverage: @{minCoverage}\n\n"))

## seems in previous step in the pipeline, chrX and chrY have beem transformed to chr23 and chr24
 
sex = read.table(gender, header=FALSE, stringsAsFactors=FALSE)[,1]
if (sex == "male" | sex == 'klinefelter') {
	chromosomes = c(1:24)
} else if (sex == "female") {
	chromosomes = c(1:23)
}

#data = read.table(file, sep = "\t", header = FALSE, as.is = TRUE, stringsAsFactors = TRUE)
colNamesData = c("chromosome", "SNP",     "start",   "end",     "crest",
                   "copyT",      "covT",    "meanTCN", "betaT",   "betaN",
                   "Atumor",     "Btumor",  "Anormal", "Bnormal", "haplotype", "map")
#dataAll = data

dataAll = lapply(chromosomes, function(chr){
		cat( "Reading chr ", chr," from ", file, "...\n" )
		dataList.chr <- try( read.table( pipe( paste( "tabix", file, chr ) ), header=FALSE, sep='\t' )  )
		if ( is.data.frame(dataList.chr)  ){
			colnames( dataList.chr ) = colNamesData
			dataList.chr
		}else{
			cat(chr," not found in ", file, "skipping!\n")
			NULL
		}
	})


segments = read.table(segments, sep = "\t", as.is = TRUE, header = TRUE, stringsAsFactors=FALSE)
segAll = segments

#get average coverage in control for maximum coverage limit of SNPs (use autosomes only)
coveragePerChr <- lapply( dataAll[1:22], function(i){ 
                                              if( ! is.null(i)){
                                                    data.frame(totalReads=sum(i[ , c('Anormal', 'Bnormal') ]), nbrSNPs=nrow(i) )
                                              }else
                                                    NULL
                                        })
totalCoverage <- do.call(rbind, coveragePerChr)
averageCoverage <- sum(as.numeric(totalCoverage$totalReads))/sum(as.numeric(totalCoverage$nbrSNPs))
cat("average SNP coverage observed in control: ", averageCoverage, "\n")

#get area under the curve for each segment's BAF distribution
segAll$area = runTheStuff(segments, chromosomes = chromosomes, minLim=minLim, maxLim=maxLim, averageCoverage, minCov=minCoverage)
cat(segAll$area[1:10],"\n")
#estimate cut off for classification  into imbalanced and balanced
densDiff <- density(segAll$area,from=0.25,to=0.35, na.rm=T)
minDiffIndex <- which(densDiff$y==min(densDiff$y))[1]
densDiffComplete <- density(segAll$area, na.rm=T)
minDiff <- which(densDiffComplete$x < densDiff$x[minDiffIndex])

png(paste0(out, "/", "peakArea_distribution.png"),width=1000,height=300, type='cairo')
 
plot( densDiffComplete, xlim=c(-1, 1), main=paste0("seperator at: ", round(densDiff$x[minDiffIndex],digits=5) ) )
abline(v=c(0.25,0.35), lty='dotted')
abline(v=densDiff$x[minDiffIndex], col='red')

dev.off()



segAll$peaks <- 1
if( any( is.na(segAll$area) ) ){
  segAll$peaks[is.na(segAll$area)] <- "missing"
}
  
selImbalanced <- which(segAll$area > densDiff$x[minDiffIndex])

if(length(selImbalanced)>1){
  segAll$peaks[selImbalanced] <- 2
}
#write.table(segAll, file = qq("@{out}/peaks.txt"), sep = "\t", row.names = FALSE, quote = FALSE) 

## it may affect, they were in `for (fac in facts)` loop
#table = subset(table, !is.na(Btumor) & !is.na(Atumor))

df.cov = vector("list", 24) 

for (chr in chromosomes) {
	#skip if no data found for chromosome and dataAll==NULL
	if ( ! is.data.frame( dataAll[[chr]] ) ){
		next
	}
	cat(qq("processing @{chr}\n\n"))
	selS = which( segAll$chromosome == chr )
	## seems unique(table$start[sel]) does not work
	facts = segAll$start[selS]
	ends  = segAll$end[selS]

	tab = dataAll[[chr]][, c("chromosome", 'SNP', "start", "covT", "Btumor", "Atumor", "betaN", "betaT"), drop = FALSE]
	tab = subset(tab, ! is.na(Btumor) & ! is.na(Atumor))
	
	df.sub = vector("list", length(facts))

	for (i in seq_along(facts)){
		sub = subset(tab, tab$start == facts[i])
		meanCovT = mean(sub$covT)
		
		i_hetero = which( (sub$betaN > 0.3) & (sub$betaN < 0.7) & ! is.na( sub$betaT ) )
		segAll$tcnNbrOfHets[selS[i]] <- length(i_hetero)
		segAll$tcnNbrOfSNPs[selS[i]] <- nrow(sub)
		#recalcluate dhMax for each segment
		dh = 2 * (abs(sub$betaT[i_hetero] - 0.5))
		if(length(i_hetero) > 0){
		  h = kde(dh, h = 0.05)                            #dh distribution for each segment
		  segAll$dhMax[selS[i]] = h$eval.points[round(mean(which(h$estimate == max(h$estimate))))] # dh maximum for each segment
		}else{
		  segAll$dhMax[selS[i]] = NA
		}

		meanCovB = median( sub[i_hetero, "Btumor" ]/(sub[i_hetero,'Btumor'] + sub[i_hetero,'Atumor']) )	#not randomly assigned alleles anymore ==> mean gives wrong value in case of LOH
		
		df.sub[[i]] = c( chr, as.numeric(facts[i]), meanCovT, meanCovB )

	}
	
	allSub = do.call(rbind, df.sub)
	df.cov[[chr]] = allSub
}


allCov = data.frame(do.call(rbind, df.cov), stringsAsFactors=FALSE)
#colnames(allCov) = c("chromosome", "start", "meanCovT", "meanCovB", "startMin", "startMax", "endMin", "endMax")
colnames(allCov) = c("chromosome", "start", "meanCovT", "meanCovB")

combi = merge(segAll, allCov, by = c("chromosome", "start"), all.x = TRUE, all.y = FALSE,)
#write.table(combi, file = paste(out, "/combi.txt", sep=''), sep="\t", row.names = FALSE, quote = FALSE)

# level covT and covB for each cluster
cluster = unique(combi$cluster[which(!is.na(combi$cluster))])
for (i in seq_along(cluster)) {

	sel = which(  combi$cluster == cluster[i] & ! is.na(combi$meanCovT) )  ## it is the same `cluster[i]` or i, but using `cluster[i]` is more safe
	selB = which( combi$cluster == cluster[i] & ! is.na(combi$meanCovB) )
	if (length(sel) < 1) {
		next
	} else {	
		combi$meanCovT[sel] = sum( as.numeric(combi$length[sel] ) * as.numeric(combi$meanCovT[sel])) / sum(as.numeric(combi$length[sel]))
#		combi$meanCovB[selB] = sum( as.numeric(combi$length[selB] ) * as.numeric(combi$meanCovB[selB])) / sum(as.numeric(combi$length[selB]))
	}
	if (length(selB) < 1) {
		next
	} else {	
		combi$meanCovB[selB] = sum( as.numeric(combi$tcnNbrOfHets[selB] ) * as.numeric(combi$meanCovB[selB])) / sum(as.numeric(combi$tcnNbrOfHets[selB]))
		combi$meanCovB[selB] = as.numeric(combi$meanCovB[selB]) * as.numeric(combi$meanCovT[selB])/10000    #divide by 10000 as coverage was initially given for 10kb windows
	}
} 

#adjust frequencies for chromosomes without cluster
NAcluster = which( is.na(combi$cluster) & ! is.na(combi$meanCovB) & ! is.na(combi$meanCovT) )
if (length(NAcluster)>0){
  combi$meanCovB[NAcluster] <- as.numeric(combi$meanCovB[NAcluster]) * as.numeric(combi$meanCovT[NAcluster])/10000
}


write.table(combi, file = segOut, sep = "\t", row.names = FALSE, quote = FALSE)
