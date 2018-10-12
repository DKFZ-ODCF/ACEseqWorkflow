#!/usr/bin/R

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

library(mclust)
library(fpc)
library(reshape)
library(ks) 
library(getopt)
library(ggplot2)
library(gridExtra)

script_dir = dirname(get_Rscript_filename())

wd = getwd()
min_num_SNPs =15
libloc=NULL
spec <- matrix(c('file',                'f', 1, "character",
                 'segments',            's', 1, "character",
				 'functions',		  	'p', 1, "character",
				 'blockPre',		  	'b', 1, "character",     # prefix of file containing haplotype groups
				 'blockSuf',		  	'u', 1, "character",     # suffix of file containing haplotype groups
				 'adjustAlleles',	  	'a', 1, "character",     # function to swap alleles where necessary
				 'sex',		  	  	  	'g', 1, "character",     # sex of patient
				 'newFile',				'w', 1, "character",
                 'out',                 'x', 1, "character",
				 'segOut',		  		'o', 1, "character",
                 'min_seg_length',      'l', 1, "numeric",
                 'clustering_YN',       'c', 1, "character",
                 'min_num_cluster',     'n', 1, "numeric",
                 'min_num_SNPs',        'i', 2, "numeric",
		 		 'min_distance',	  	'd', 1, "numeric",
                 'min_membership',      'm', 1, "numeric",
				 'chrLengthFile',       'r', 1, "character",
				 'gcCovWidthFile',      'y', 1, "character",
				 'pid',			  		'v', 1, "character",
				 'libloc',		  		'z', 2, 'character',
				 'runInDebugMode',		'k', 2, 'character'
                ), ncol = 4, byrow = TRUE)

opt = getopt(spec);
for(item in names(opt)){
       assign( item, opt[[item]])
}


cat(paste0("file: ",file, "\n\n"))
cat(paste0("segments: ",segments, "\n\n"))
cat(paste0("functions: ",functions, "\n\n" ))
cat(paste0("out: ",out, "\n\n"))
cat(paste0("segOut: ",segOut, "\n\n"))
cat(paste0("blocks: ",blockPre, "*",blockSuf, "\n\n"))
cat(paste0("adjustAlleles: ",adjustAlleles, "\n\n"))
cat(paste0("gcCovWidthFile: ",gcCovWidthFile, "\n\n"))
cat(paste0("newFile: ",newFile, "\n\n"))
cat(paste0("min_seg_length: ",min_seg_length, "\n\n"))
cat(paste0("clustering_YN: ",clustering_YN, "\n\n"))
cat(paste0("min_num_cluster: ",min_num_cluster, "\n\n"))
cat(paste0("min_num_SNPs: ",min_num_SNPs, "\n\n"))
cat(paste0("min_distance: ",min_distance, "\n\n"))
cat(paste0("min_membership: ",min_membership, "\n\n"))
cat(paste0("chrLengthFile: ",chrLengthFile, "\n\n"))
cat(paste0("gcCovWidthFile: ",gcCovWidthFile, "\n\n"))
cat(paste0("runInDebugMode: ",runInDebugMode, "\n\n"))
cat("\n")


if ( libloc == "" | libloc == TRUE ){
	libloc=NULL
}

source(functions)

mclust.options(hcUse = "VARS") # set hcUse to use VARS-method (method will be changed to SVD mith mclust>=5.4 which leads to different results)

cat(paste0("reading ",segments,"...\n\n"))
segAll = read.table(segments, sep = "\t", as.is = TRUE, header = TRUE)

#########################
### remove umapped
segAll = segAll[segAll$map != "unmappable", ]

#read in limits of main Cluster
#fieldsIn_FILENAME_GC_CORRECTED_QUALITY=c(pid, half_max_pos_left, half_max_pos_right, GCcorrectQuant_string, mean_normal_slope, mean_abs_normal_slope,mean_tumor_slope,
#mean_abs_tumor_slope, mean_normal_curvature, mean_abs_normal_curvature,mean_tumor_curvature, mean_abs_tumor_curvature,
#main_cluster_width_n, main_cluster_width_t, main_cluster_FWHM_n, main_cluster_FWHM_t, mean_abs_delta_slope,
#mean_delta_slope, mean_abs_delta_curvature, mean_delta_curvature, minimal_coverage_gcfit_normal, minimal_coverage_gcfit_tumor ),
#sep="\t",file=outGCcorrectQuant_file, ncolumns=37 )
covWidthLimits <- read.table(gcCovWidthFile, header=F)[,2:3] #half_max_pos_left, half_max_pos_right
covLeft <- covWidthLimits[,1] #half_max_pos_left
covRight <- covWidthLimits[,2] #half_max_pos_right
covWidth <- covRight -covLeft

sex <- read.table(sex, header=FALSE, stringsAsFactors=FALSE)[,1]
cat( paste0("sex: ", sex, "\n\n") )

if (sex =='male'){
	#in case the patient is male do not take X and Y chromosome into account as they are haploid
	maxChr = 24
	chromosomes = c(1:22)
	segXY  <- segAll[ segAll$chromosome==23 | segAll$chromosome==24, ]
	segAll <- segAll[ segAll$chromosome!=23 & segAll$chromosome!=24, ]
}else if (sex== 'klinefelter'){
	maxChr = 24
	chromosomes = c(1:23)
	segXY  <- segAll[ segAll$chromosome==24, ]
	segAll <- segAll[ segAll$chromosome!=24, ]
}else {
	maxChr = 23
	chromosomes = c(1:23)
}

cat(paste0("blocks: ",blockPre,"*",blockSuf,"\n\n"))
for (chr in chromosomes) {
	blockFile <- NULL
	blockFile <- paste0( blockPre, chr, ".", blockSuf)
	fileCheck <- file.exists(blockFile)
	if ( fileCheck == FALSE ) {
		cat( "haploblock file for ", chr, " not found!! Exiting..." )
		quit( save="no", status=2 )	
	}
}

#read chromosome length file
chr = read.table(chrLengthFile, header = FALSE, as.is = TRUE)
chrLength = data.frame(chr)
chrLength$V1 <- gsub('chr','',chrLength$V1)
sel = which(chrLength$V1 == "X")
chrLength$V1[sel] = 23
sel = which(chrLength$V1 == "Y")
chrLength$V1[sel] = 24


colNamesData <-  c( "chromosome", "SNP", "start", "end", "SV.Type", "copyT", "covT", "meanTCN", "betaT","betaN", "Atumor", "Btumor", "Anormal", "Bnormal", 'haplotype', "map" )

dataAll = lapply( seq_len(maxChr), function(chr){
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



# delete, attach or ignore very short segments
if (min_seg_length != 0) {
	
	cat("delete, attach or ignore very short segments\n\n")
	rows = rep(1:nrow(segAll))
	nrows = nrow(segAll)
	i = 0
	for (r in seq_len(nrows)) {
		
		i = i + 1
		#first segment in chromosome
		if ( i != nrows && ( (segAll$length[i] < min_seg_length ||  segAll$tcnNbrOfHets[i] < min_num_SNPs ) && segAll$map[i] != 'homozygousDel'  ) && 
		    (i == 1 || segAll$chromosome[i] != segAll$chromosome[i-1]) &&
			 ( segAll$end[i]==segAll$start[i+1] | segAll$end[i]+1==segAll$start[i+1] ) ) { 
		        
			if ( segAll$chromosome[i] != segAll$chromosome[i+1] | segAll$map[i+1]=="homozygousDel"  )  { 

				cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment is smaller min_seg_length and only one in chromosome -> not attached to next segment\n\n"))
				next
			}

			cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment is smaller min_seg_length and first in chromosome -> attached to next segment\n\n"))
			
			segAll$start[i + 1] = segAll$start[i]
			segAll$length[i + 1] = segAll$length[i + 1] + segAll$length[i]
			segAll = segAll[-i, , drop = FALSE]

			i = i-1
			nrows = nrows-1
		
    #last segment in chromosome
		} else if ( i != 1 && ( (segAll$length[i] < min_seg_length ||  segAll$tcnNbrOfHets[i] < min_num_SNPs ) && segAll$map[i] != 'homozygousDel' ) &&
			    ( i == nrows || segAll$chromosome[i] != segAll$chromosome[i+1] ) &&
				 ( segAll$end[i-1] == segAll$start[i] | segAll$end[i-1]+1==segAll$start[i] ) ) {

			if (  segAll$chromosome[i] != segAll$chromosome[i-1] | segAll$map[i-1]== "homozygousDel" ) {

				cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment is smaller min_seg_length and only one in chromosome -> not attached to next segment\n\n"))
				next
			}
				
			cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment is smaller min_seg_length and last in chromosome -> attached to former segment\n\n"))
			
			segAll$end[i - 1] = segAll$end[i]
			segAll$length[i - 1] = segAll$length[i - 1] + segAll$length[i]
			segAll = segAll[-i, , drop = FALSE]
			
			i = i-1
			nrows =  nrows-1
			
		} else if (i != 1 && i != nrows && 
		           ( (segAll$length[i] < min_seg_length ||  segAll$tcnNbrOfHets[i] < min_num_SNPs ) && segAll$map[i] != 'homozygousDel' )  && 
				   segAll$chromosome[i] == segAll$chromosome[i - 1] && 
				   segAll$chromosome[i] == segAll$chromosome[i + 1] ) {
		
			cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment is smaller min_seg_length and within chromosome\n\n"))
			
      #segment is only adjacent to prior segment or more similar to prior segment with regards to coverage
			if (i != 1 && i != nrows && 
			    is.na(segAll$SV.Type[i]) &&
				 ( segAll$end[i-1] == segAll$start[i] | segAll$end[i-1]+1==segAll$start[i] ) &&
           segAll$map[i-1]!="homozygousDel" &&
					( (segAll$start[i+1] != segAll$end[i] & segAll$end[i]+1 != segAll$start[i+1] ) ||
					  abs(segAll$tcnMean[i] - segAll$tcnMean[i - 1]) <= abs(segAll$tcnMean[i] - segAll$tcnMean[i + 1]) ) ) {
			
				cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment start not eq sv point and tcnMean differnence to former segment < to next segment -> attached to former segment\n\n"))
			
				segAll$end[i - 1] = segAll$end[i]
				segAll$length[i - 1] = segAll$length[i - 1] + segAll$length[i]
				segAll = segAll[-i, , drop = FALSE]
				
				i = i - 1
				nrows = nrows - 1
			
      #segment is only adjacent to following segment or more similar to following segment with regards to coverage  
			} else if (i != 1 && i != nrows && 
			           is.na(segAll$SV.Type[i+1]) && 
                 segAll$map[i+1]!="homozygousDel" &&
					( segAll$end[i] == segAll$start[i+1] | segAll$end[i]+1==segAll$start[i+1] ) &&
						( (segAll$start[i] != segAll$end[i-1] & segAll$end[i-1]+1 != segAll$start[i] ) ||
						   abs(segAll$tcnMean[i] - segAll$tcnMean[i - 1]) > abs(segAll$tcnMean[i] - segAll$tcnMean[i + 1]) ) ) {
			
				cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment start not eq sv point and tcnMean differnence to former segment > to next segment -> attached to next segment\n\n"))
			
				segAll$start[i + 1] = segAll$start[i]
				segAll$length[i + 1] = segAll$length[i + 1] + segAll$length[i]
				segAll = segAll[-i, , drop = FALSE]
				
				i = i - 1
				nrows = nrows - 1
				
			} else {
			
				cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment has sv point -> not attached\n\n"))
				next
				
			} 
			
		} else {
		
			cat(paste0("",segAll$chromosome[i], "-",segAll$start[i], ": segment longer min_seg_length -> not attached\n\n"))
			next
			
		}
    }
    
}
       
# determine maximum of dh distribution for each segment
cat("determine maximum of dh distribution for each segment\n\n")

                         
i = 0
dhMax = seq_len(nrow(segAll)) 
                         
for (chr in chromosomes) {
	
    cat(paste0("processing ",chr, "...\n\n"))
    
    dataAll[[chr]]$dh = 2 * (abs(dataAll[[chr]]$betaT - 0.5)) # dh value for each SNP
    sel = which(segAll$chromosome == chr)
    start = segAll$start[sel]		#the snps must be chosen according to SNP location not start, as it could be that relevant SNPs are missed otherwise!
    end = segAll$end[sel]
    j = 0
    
    for (seg in start) {
	
	    i = i + 1
	    j = j + 1
			
	    selLoci <- which( dataAll[[chr]]$SNP >= start[j] &
			      dataAll[[chr]]$SNP <= end[j] )
	    
	    segAll$tcnNbrOfLoci[i] <-  length(selLoci) 
	      
	    s = which(  dataAll[[chr]]$betaN[selLoci] > 0.3 & 
			dataAll[[chr]]$betaN[selLoci] < 0.7 & 
                      ! is.na(dataAll[[chr]]$dh[selLoci]) )
			
	    if (length(s) > 0) {
	      h = kde(dataAll[[chr]]$dh[selLoci[s]], h = 0.05)	#dh distribution for each segment
	      dhMax[i] = h$eval.points[round(mean(which(h$estimate == max(h$estimate))))] # dh maximum for each segment
	      segAll$tcnNbrOfHets[i] <- length(s)
	    } else {
	      dhMax[i] = NA
	      segAll$tcnNbrOfHets[i] <- 0
	    }     
  }
}
                  
dhMax_new = as.numeric(dhMax)
segAll$dhMax = dhMax_new
write.table(segAll, file = paste0("",out, "/minus_short_w_dh_BIC.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

tcnMean = segAll$tcnMean
dhMax = segAll$dhMax

png(paste0("",out, "/",pid, "_tcn_dh.png"), width=1200, height=1200, type='cairo')
plot(tcnMean,dhMax)
dev.off()
  
if (clustering_YN == "yes") {
	cat("clustering YN...\n")
	# cluster_matrix is matrix of normalized/scaled tcnMean and dhMax values 
	cluster_matrix_norm = cbind(log2.tcnMean=log2(tcnMean),dhMax= dhMax)

	rem = which(is.na(cluster_matrix_norm[,1]) | 
	            is.na(cluster_matrix_norm[,2]) | 
	            segAll$map == "unmappable" | 
	            cluster_matrix_norm[ ,1] == "Inf" | 
	            cluster_matrix_norm[ ,1] == "-Inf" | 
	            cluster_matrix_norm[ ,2] == "Inf" | 
	            cluster_matrix_norm[,2] == "-Inf")

	cat(paste0("",length(rem), " lines dropped.\n\n")) #
	if ( length(rem) > 0) {
		# +1 as pseudo count to account for segments of length 1bp
		weights = log2(segAll$length[-rem] + 1) #specify your size vector here
												 #weights <- 1
		cluster_matrix_norm = cluster_matrix_norm[-rem, ]
		#convert limits to scaled coordinates
		covLeftNorm  <- ( log2(covLeft)-mean(log2(tcnMean[-rem])) )/sd(log2(tcnMean[-rem]))
		covRightNorm <- ( log2(covRight)-mean(log2(tcnMean[-rem])) )/sd(log2(tcnMean[-rem]))
		covLeftFullNorm <- ( log2(covLeft-covWidth/2)-mean(log2(tcnMean[-rem])) )/sd(log2(tcnMean[-rem]))
		covRightFullNorm <- ( log2(covRight+covWidth/2)-mean(log2(tcnMean[-rem])) )/sd(log2(tcnMean[-rem]))

  } else {
		# +1 as pseudo count to account for segments of length 1bp
		weights = log2(segAll$length + 1)
		cluster_matrix_norm = cluster_matrix_norm
		#convert limits to scaled coordinates
		covLeftNorm  <- (log2(covLeft)-mean(log2(tcnMean)) ) / sd(log2(tcnMean))
		covRightNorm <- (log2(covRight)-mean(log2(tcnMean)) ) / sd(log2(tcnMean))
		covLeftFullNorm <- ( log2(covLeft-covWidth/2)-mean(log2(tcnMean)) ) / sd(log2(tcnMean))
		covRightFullNorm <- ( log2(covRight+covWidth/2)-mean(log2(tcnMean)) ) / sd(log2(tcnMean))
  }
  
	cluster_matrix = scale(cluster_matrix_norm)

	#find optimal number of clusters using bayesian information criterion
	cat(paste0(Sys.time(),": Calling Mclust...\n"))
	cat(paste0("cluster_matrix nrow: ",nrow(cluster_matrix),"\n"))
	cat(paste0("cluster_matrix ncol: ",ncol(cluster_matrix),"\n"))
	cat(paste0("min_num_cluster: ",min_num_cluster,"\n"))
	d_clust <- Mclust(cluster_matrix, G=min_num_cluster:20)
	cat(paste0(Sys.time(),": finished Mclust...\n"))
	m.best  <- dim(d_clust$z)[2]

	#cmeans to get clusters with m.best centers 
	#resample clustering by jittering point B times)

	cat(paste0(Sys.time(),": Calling clusterboot...\n"))
	results = clusterboot(cbind(weights, cluster_matrix), B = 100, bootmethod = "jitter", clustermethod = cmeansCBI, k = m.best, seed = 15555, multipleboot = FALSE)
	cat(paste0(Sys.time(),": finished clusterboot...\n"))

	if ( runInDebugMode == "true") {
		save.image(paste0("",out, "/",pid, "_cluster_data.RData"))
	}

	CM <- results$result$result

    	massCenterX <- sapply(seq_along(CM$centers[,1]), function(i) {
    	  median(cluster_matrix[CM$cluster==i,1])
    	})
    
    	massCenterY <- sapply(seq_along(CM$centers[,1]), function(i) {
      		median(cluster_matrix[CM$cluster==i,2])
    	})
    
    	massCenter <- cbind(massCenterX, massCenterY)

	frequencies <- table(CM$cluster)
	clusterWithinLimits <- which(CM$centers[,1] > covLeftNorm & CM$centers[,1] < covRightNorm )
  centerMain <- NULL
  if (length(clusterWithinLimits) >1 ){
  #in case several samples fullfill criteria with equal amount of points
    set.seed(seed=15555)
    maxCluster <- sample( names( which(frequencies[clusterWithinLimits]==max(frequencies[clusterWithinLimits]) ) ), size= 1 )
	  maxCluster <- as.numeric(maxCluster)
	  centerMain <- CM$centers[maxCluster,]
  }
	col = c("#000000","#800000","#008000","#000080","#800080","#808080","#FF0000","#00FF00","#FFFF00","#0000FF","#FF00FF","#00FFFF","#DC143C","#FF8C00","#FF69B4","#FF4500", "#EE82EE", "#FFD700")

	clusterPlot <- ggplot( data.frame(cluster_matrix), aes(log2.tcnMean, dhMax, col=as.character(CM$cluster) ) )  +
	  geom_point(size=1.7, alpha=0.8) +
	  geom_point(data=data.frame(CM[['centers']]), aes(log2.tcnMean, dhMax), col='red', pch=3) +
	  geom_vline(xintercept=c(covRightNorm, covLeftNorm),size=0.4, col="black", alpha=0.8) +
	  geom_vline(xintercept=c(covRightFullNorm, covLeftFullNorm),size=0.4, col="red", alpha=0.8) +
	  scale_color_manual(values=c(col[1:length(unique(CM$cluster))], "grey"), name="cluster" )
	ggplot2::ggsave(paste0("",out, "/",pid, "_cluster_cmeans.png"), clusterPlot, width = 10, height = 10, type='cairo')

	minTcnMean <- covLeftNorm
	maxTcnMean <- covRightNorm
	
  if ( ! is.null(centerMain)){
	  CM_merged <- mergeClusters(CM, minTcnMean, maxTcnMean,cluster_matrix, maxCluster)
  }else{
    CM_merged <- CM
  }
	
  massCenterX <- sapply(seq_along(CM_merged$centers[,1]), function(i) {
	  median(cluster_matrix[CM_merged$cluster==i,1])
	})
	
	massCenterY <- sapply(seq_along(CM_merged$centers[,1]), function(i) {
	  median(cluster_matrix[CM_merged$cluster==i,2])
	})
	
	massCenter <- cbind(massCenterX, massCenterY)
  
	CM_new <- CM_merged
  
  segAll$cluster = rep(0, nrow(segAll))
  for (i in seq_along(rem)) {
    segAll$cluster[rem[i]] <- NA
  }
  
  sel = which(!is.na(segAll$cluster))
  for (i in seq_along(CM_new$cluster)) {
    segAll$cluster[sel[i]] = CM_new$cluster[i]
  }     
  
  segAll.tmp <- combineNeighbours(segAll)
  keep <- (! is.na(segAll.tmp$tcnMean) & ! is.na(segAll.tmp$dhMax))
  
  #set clusters with less than 5 members to NA
  frequencies <- table(segAll.tmp$cluster)
  minimalClusters <- names(which(frequencies < 5 ))
  for(i in as.numeric(minimalClusters)){
    selMin <- which(segAll.tmp$cluster==i )
    segAll.tmp$cluster[selMin] <- NA
  }
  
  #determine centers by the point of highest density in x and y direction
  newCenters <- data.frame( t(sapply(1:max(unique(CM_new$cluster)), function(i){
    s <-  which(segAll.tmp$cluster==i)
    tcnMean <- NaN
    dhMax <- NaN
    if(length(s)>4){
      tcnMean <- density(segAll.tmp$tcnMean[s] )
      tcnMean <- tcnMean$x[which(tcnMean$y==max(tcnMean$y))[1]] # [1]: bugfix, more than one value possible. always take the first one
      dhMax <- density(segAll.tmp$dhMax[s])
      dhMax <- dhMax$x[which(dhMax$y==max(dhMax$y))[1]] # [1]: bugfix, more than one value possible. always take the first one
    }
    c(tcnMean, dhMax)
  })))
  
  colnames(newCenters) <- c("tcnMean", "dhMax")
  segAll.tmp[keep,] <- removeOutlierPoints_cmean_alt( segAll.tmp[keep,],  newCenters, deviationFactor = 2)
  
  clusterPlotRmOutlier  <- ggplot( data.frame(segAll.tmp), aes(log2(tcnMean), dhMax, col=as.factor(cluster) ) )  +
    geom_point(size=1.7, alpha=0.8) +
    geom_vline(xintercept=c(log2( covLeft) , log2(covRight) ), size=0.4, col="black", alpha=0.8) +
    geom_vline(xintercept=c(log2( covLeft -  covWidth), log2(covRight + covWidth) ),size=0.4, col="red", alpha=0.8) +
    scale_color_manual(values=c(col[1:(length(unique(segAll.tmp$cluster))-1)]),name="cluster" ) +
    geom_point(data=data.frame(segAll.tmp[is.na(segAll.tmp$cluster),]), aes(log2(tcnMean), dhMax), col='grey')
	ggplot2::ggsave(paste0("",out, "/",pid, "_cluster_cmeans_wo_outlier.png"), clusterPlotRmOutlier, width=10, height=10, type='cairo')
  
  segAll.tmp[keep,]$cluster[is.na(segAll.tmp[keep,]$cluster)] <- "NA"  #seem redundant and uneccesary
 
  freqs <- sapply( as.numeric(rownames(newCenters)),
               function(i) sum(segAll.tmp$cluster==i, na.rm=T) )
  names(freqs) <- rownames(newCenters)                     
  minimalClusters <- names(which(freqs <5 ))
  for(i in as.numeric(minimalClusters)){
	  selMin <- which(segAll.tmp$cluster==i )
	  segAll.tmp$cluster[selMin] <- NA
          newCenters[as.character(i),] <- NA
  }
  # write.table(segAll, file = paste0("",out, "/clustered.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # pruning if neighbouring segments are assigned to identical cluster                
	
  test <- combineNeighbours(segAll.tmp)

  if ( any( is.na(newCenters)) ) {
          newCenters <- newCenters[! is.na(newCenters[,1]),]
  }
  #name x direction for mainCluster estimation
  newCentersTcnMean <- newCenters[,1]
  names(newCentersTcnMean) <- rownames(newCenters)

  frequencies <- table(test$cluster)
  clusterWithinLimits <- names(which( newCentersTcnMean > covLeft & newCentersTcnMean < covRight ))
  centerMain <- NULL
  if (length(clusterWithinLimits) >0 ){
    set.seed(seed=15555)
    maxCluster <- sample( names( which(frequencies[as.character(clusterWithinLimits)]==max(frequencies[as.character(clusterWithinLimits)])) ), size= 1 )
    maxCluster <- as.numeric(maxCluster)
    centerMainNew <- data.frame( newCenters[as.character(maxCluster),])

    test$neighbour <- findNeighbours(test, maxCluster)
    selIdentical <- which(test$neighbour=='identical')
    
    test_new <- mergePoints(test,centerMainNew,maxCluster, covLeft, covRight)
    test_new$neighbour <- findNeighbours(test_new, mainCluster=maxCluster)
  
  }else{
    test_new <- test
    test_new$neighbour <- NA
    test_new$distDH <- NA
    test_new$errorSNP <- NA
    test_new$distTcn <- NA
    test_new$errorLength <- NA
    test_new$totalError <- NA
  }
  
  #revert outlier to real NAs
  if (any (test_new$cluster == "NA",na.rm=T) )
	  test_new$cluster[test_new$cluster=="NA"] <- NA
    selIdentical <- which( test_new$neighbour == 'identical' )

    minNbrOfHets =5
  png(paste0("",out, "/",pid, "_merged_cluster.png"), width=1500, height=1020, type='cairo')
    generatePlots(test_new, selIdentical, centerMainNew, covLeft, covRight, covLeft-covWidth, covRight+covWidth, minNbrOfHets=minNbrOfHets)
  dev.off()
      

  clusterPlotNewlog2 <- ggplot( data.frame(test_new), aes(log2(tcnMean), dhMax, col=as.factor(cluster) ) )  +
    geom_point(size=1.7, alpha=0.8) +
    geom_point(data=data.frame(test_new[is.na(test_new$cluster),]), aes(log2(tcnMean), dhMax), col='grey') +
    geom_vline(xintercept=c(log2( covLeft) , log2(covRight) ), size=0.4, col="black", alpha=0.8) +
    geom_vline(xintercept=c(log2( covLeft -  covWidth), log2(covRight + covWidth) ),size=0.4, col="red", alpha=0.8) +
    scale_color_manual(values=c(col[1:(length(unique(test_new$cluster))-1)]),name="cluster" )
	ggplot2::ggsave( paste0("",out, "/",pid, "_cluster_cmeans_merged_log2.png"), clusterPlotNewlog2, width=10, height=10, type='cairo' )
#	write.table(test, file = paste0("",out,"/clustered_and_pruned_BIC.txt"), sep = "\t", row.names = FALSE, quote = FALSE )

  combi <- test_new

  chrCount <- max(combi$chromosome)
  chrNames <- as.character(1:chrCount)
  maxCov = max(combi$tcnMean, na.rm=TRUE)
  
  #plot per chromosomes
  for (chr in seq_len(chrCount)){
    selSeg <- which(combi$chromosome==chr)
    if (length(selSeg) > 0){
      chrL = as.numeric(chrLength[chrLength$V1==chr,2]) 
      pChr <- plotCov(combi[selSeg,],chrLen=chrL)
      pChr1 <- pChr[[1]]
      pChr2 <- pChr[[2]]
      hlinesDH <- geom_hline( yintercept= c(0, 0.2, 0.4, 0.6, 0.8, 1), col="lightgray", lty="dotted", lwd=0.4)
      pChr1 <- pChr1 + geom_hline( yintercept=seq(0, maxCov+0.8, 1), col='lightgray', lty='dotted', lwd=0.4 )
      pChr1 <- pChr1 + labs(x=NULL) + theme( axis.text.x=element_blank() )
      
      pChr2 <- pChr2 + hlinesDH
      pChr2 <- pChr2 + labs(x=NULL) + theme( axis.text.x=element_blank() )
      pC <- arrangeGrob(pChr1,pChr2)
      ggplot2::ggsave( paste0(out, "/", pid,"_chr", chr,".pdf"), pC, width=28,height=14 )
    }  
  }
  
  chrL = sum( as.numeric(chrLength[1:chrCount,2]) )
  
  xoff  = 0 
  xoffsets = xoff
  xtot = chrL/10
  
  #adjust coordinates so plots can be made from single dataframe with single command
  for ( chr in seq_len(chrCount) ) { 
    
    sel = which(combi$chromosome==chr)
    combi$start[sel] <- combi$start[sel]+xoff
    combi$end[sel]   <- combi$end[sel]+xoff
    
    xoff = xoff + chrLength[chrLength$V1==chr,2] 
    xoffsets = append( xoffsets, xoff/10 )
  }   
  
  p  <-  plotCov(combi,chrLen=chrL)
  p1 <- p[[1]]    
  p2 <- p[[2]]
  
  #add chromosome labels and borders
  labels.data        <- data.frame(xPos=chrLength[1:chrCount,2] / 20 + xoffsets[1:length(xoffsets) -1] )
  labels.data$xPos   <- labels.data$xPos/xtot
  labels.data$maxCov <- replicate(nrow(labels.data), maxCov)
  labels.data$chr    <- as.character(chrNames)
  
  hlinesDH <- geom_hline( yintercept= c(0, 0.2, 0.4, 0.6, 0.8, 1), col="lightgray", lty="dotted", lwd=0.4)
  vlines   <- geom_vline( xintercept= xoffsets[2:length(xoffsets)]/xtot, col="#000000", lwd = 0.2 )
  
  p1 <- p1 + geom_hline( yintercept=seq(0, 2, 0.5), col='lightgray', lty='dotted', lwd=0.4 )
  p1 <- p1 + vlines + labs(x=NULL) + theme( axis.text.x=element_blank() )
  p1 <- p1 + geom_text(data=labels.data, aes(x=xPos, y=1.4, label=chr), size=8)
  p2 <- p2 + hlinesDH
  p2 <- p2 + vlines + labs(x=NULL) + theme( axis.text.x=element_blank() )
  p <- arrangeGrob(p1,p2)
  ggplot2::ggsave(paste0(out, "/", pid, "_colored_segments_cluster.pdf"), p,width=28,height=14,units='cm' ) 
 # }


	# same tcnMean for each segment of one cluster          
	cluster_tcn = unique(test_new$cluster)[which(!is.na(unique(test_new$cluster)))]
	tcnMean_sum = rep(0, length(cluster_tcn))
	for (i in seq_along(cluster_tcn)) {
	
		sel = which(test_new$cluster == cluster_tcn[i])
		tcnMean_new = rep(0, length(sel))
		loci_new = rep(0, length(sel))
		
		for (j in seq_along(sel)) {
			tcnMean_new[j] = test_new$tcnNbrOfLoci[sel[j]]*test_new$tcnMean[sel[j]]
			loci_new[j] = test_new$tcnNbrOfLoci[sel[j]]
		}
		
		tcnMean_sum[i] = sum(tcnMean_new) / sum(loci_new)
		test_new$tcnMean[sel] = tcnMean_sum[i]
	}

} else if (clustering_YN == "no") {
	segAll$cluster = "NA"
	test_new = segAll
}


if ( sex == 'male' | sex == 'klinefelter' ){
	if ( nrow(segXY) > 0 ){

		segXY$dhMax <- NA
		segXY$cluster <- NA
		segXY$distDH <- NA
		segXY$errorSNP <- NA
		segXY$distTcn <- NA
		segXY$errorLength  <- NA
		segXY$totalError <- NA
		segXY$neighbour <- NA

		test_new <- rbind( test_new, segXY )
	}
}

test_new$minStart =NA 
test_new$maxStart =NA 
test_new$minStop =NA 
test_new$maxStop =NA
    
#add allele adjustment here (use a function that returns dataAll?
source(adjustAlleles)
 
for (chr in  seq_len(maxChr) ) {

 	cat (paste0("adjusting frequencies for chromosome ",chr, " \n\n") )
 	selSeg <- which(test_new$chromosome==chr)
 
 	if ( length(selSeg) > 0  & is.data.frame(dataAll[[chr]]) ){
		#find first and last SNP of previous and following segment
		starts	= unique(test_new$start[selSeg])	 #start points of all segments in chr
		ends 	  = unique(test_new$end[selSeg])
#		starts	= sort(unique(starts))

		for (i in seq_along(starts)){

			subSNP <- which(dataAll[[chr]]$SNP >= starts[i] & dataAll[[chr]]$SNP <= ends[i])
			if (length(subSNP) > 0 ) {
				test_new$maxStart[selSeg[i]] = suppressWarnings(min(dataAll[[chr]]$SNP[subSNP]))
				test_new$minStop[selSeg[i]]  = suppressWarnings(max(dataAll[[chr]]$SNP[subSNP]))
			}else{
				test_new$maxStart[selSeg[i]] = starts[i] 
				test_new$minStop[selSeg[i]]  = ends[i]
			}
			selPre = NULL
			selPost = NULL
			if(length(starts)>1){
				if ( i == 1 ){							#first segment on chromosome
					test_new$minStart[selSeg[i]]  = 0
					if(ends[i] == starts[i+1] | ends[i] == starts[i+1]+1){
						selPost =  which(dataAll[[chr]]$SNP >= starts[i+1] & dataAll[[chr]]$SNP <= ends[i+1])
						test_new$maxStop[selSeg[i]] = suppressWarnings( min( as.numeric(dataAll[[chr]]$SNP[selPost]) ) )
					}else
						test_new$maxStop[selSeg[i]] = ends[i]
				}else if ( i == length(starts) ){				#last segment on chromosome
					test_new$maxStop[selSeg[i]]	= ends[i]
					if(ends[i-1] == starts[i] | ends[i-1] == starts[i]+1){
						selPre =  which(dataAll[[chr]]$SNP >= starts[i-1] & dataAll[[chr]]$SNP <= ends[i-1])
						test_new$minStart[selSeg[i]]	= suppressWarnings( max( as.numeric(dataAll[[chr]]$SNP[selPre]) ) )
					}else
						test_new$minStart[selSeg[i]] = starts[i]
				
				}else{								#middle segments on chromosome
					if(ends[i-1] == starts[i] | ends[i-1] == starts[i]+1){
						selPre		=  which(dataAll[[chr]]$SNP >= starts[i-1] & dataAll[[chr]]$SNP <= ends[i-1])
						test_new$minStart[selSeg[i]]	= suppressWarnings( max( as.numeric(dataAll[[chr]]$SNP[selPre]) ) )
					}else
						test_new$minStart[selSeg[i]] = starts[i]
					if(ends[i] == starts[i+1] | ends[i] == starts[i+1]+1){
						selPost		=  which(dataAll[[chr]]$SNP >= starts[i+1] & dataAll[[chr]]$SNP <= ends[i+1])
						test_new$maxStop[selSeg[i]]	=  suppressWarnings( min( as.numeric(dataAll[[chr]]$SNP[selPost]) ) )
					}else
						test_new$maxStop[selSeg[i]] = ends[i]
				}
			}else{
				test_new$maxStop[selSeg[i]]  = ends[i]
				test_new$minStart[selSeg[i]] = starts[i]
			}

		   	 #for segments that neighbor homozygous deletions or are homozygous deletion regions
			if ( length(selPre) < 1 || test_new$minStart[selSeg[i]]==Inf || test_new$minStart[selSeg[i]] == -Inf ){
				test_new$minStart[selSeg[i]] = starts[i]
			}
			if ( length(selPost) < 1 || test_new$maxStop[selSeg[i]]==Inf || test_new$maxStop[selSeg[i]] == -Inf ){
				test_new$maxStop[selSeg[i]] = ends[i]
			}
			rm(selPre)
			rm(selPost)
		}

		if ( chr <= max(chromosomes) ) {
		#adjust allele frequencies
			dataAll[[chr]]$adjusted <- NA
			dataAll[[chr]] <- swapAlleles( test_new[selSeg,], dataAll[[chr]], chr, blockPre, blockSuf)
   
			selRem <- which( dataAll[[chr]]$betaN > 0.3 & dataAll[[chr]]$betaN < 0.7 & is.na(dataAll[[chr]]$adjusted) )
			if ( length(selRem) > 0 ){
				cat( paste0("Removing ",length(selRem), " unphased SNPs from chromosome ",chr, "...\n\n") )
				dataAll[[chr]] <- dataAll[[chr]][-selRem,]
			}
 
				dataAll[[chr]] <- dataAll[[chr]][,! colnames(dataAll[[chr]]) == 'adjusted']
				dataAll[[chr]] <- dataAll[[chr]][,! colnames(dataAll[[chr]]) == 'dh']

		}
 
   	}
}
 

 
# for chromosome 1 to change columnnames for saving file so segments to data work
new_colnames <- c("chromosome", "x", "a", "end", "SV.Type", "CT", "covT", "meanTCN", "betaT","betaN", "Atumor", "Btumor", "Anormal", "Bnormal", 'haplotype', "map")

dataAll[[1]] = format(dataAll[[1]], scientific = FALSE, trim = TRUE)
write.table(dataAll[[1]], pipe(paste0("bgzip >", newFile) ), sep='\t', col.names=new_colnames, row.names=FALSE, quote=FALSE) 

# append chr 2 - 24
writetable <- function(data, newFile){
  format = format(data, scientific = FALSE, trim = TRUE)
  write.table(format, pipe( paste0( "bgzip >>",newFile ) ), append=TRUE, sep='\t', col.names=FALSE,  row.names=FALSE, quote=FALSE)
}

lapply(dataAll[-1], writetable, newFile = newFile)


#dataAll <- do.call(rbind, dataAll)
# 
##change columnnames for saving file so segments to data works
#colnames(dataAll) = c("chromosome", "x", "a", "end", "sv", "CT", "covT", "meanTCN", "betaT","betaN", "Atumor", "Btumor", "Anormal", "Bnormal", 'haplotype', "map")
#dataAll = format(dataAll, scientific = FALSE, trim = TRUE)
#write.table(dataAll, pipe(paste0("bgzip >",newFile) ), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE) 

#write segments
test_new2 = format(test_new, scientific = FALSE, trim = TRUE)
write.table(test_new2, file = segOut, sep = "\t", row.names = FALSE, quote = FALSE )
