#!/usr/bin/R

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).

library(e1071)
library(flexclust, lib.loc=libloc)
library(reshape)
library(plyr)

# length weighted cmeans clustering with jitter noise, result: partial memberships for each segment  
cmeansCBI = function (data, krange, k = NULL, scaling = FALSE, runs = 1, criterion = "ch",...) {
  if (!is.null(k)) {
    krange <- k
  }
  
  weights = data[, 1]
  data = data[, 2:ncol(data)]
  sdata = scale(data, scale = scaling, center = scaling)
  c1 = cmeans(sdata, centers = k, weights = weights)
  partition = c1$cluster
  cl = lapply(seq_len(krange), function(i) partition == i)
  list(result = c1, nc = krange, clusterlist = cl, partition = partition, clustermethod = "cmeans")
}

findMergeBIC <- function(cluster, data, s, classes, min_distX, min_distY=0){
  
  merged <- matrix(NA, nrow=length(sel), ncol=2)
  
  #find smallest distance for each class pair
  for (r in seq_along(s)){
    clust1 <- data[cluster$cluster==classes$c1[s[r]] & ! is.na (cluster$cluster),]
    clust2 <- data[cluster$cluster==classes$c2[s[r]] & ! is.na (cluster$cluster),]
    
    # get cov R distance and merge if smaller min_dist
    if (nrow(clust1) > 3 & nrow(clust2) > 3){
      dists   <- findMinimum( dist2(clust1, clust2, method='manhattan') )
      distsX  <- dist2(clust1[,1], clust2[,1])[dists[,1:2]]
      distsY  <- dist2(clust1[,2], clust2[,2])[dists[,1:2]]
      if (sum(distsX < min_distX) == 3 & sum(distsY < min_distY) == 3  ){
        merged[r,] <- c( max(classes[s[r], 1:2]), min(classes[s[r], 1:2]) ) 
      }
    }
  }
  
  merged <- matrix(merged[ ! is.na(merged) ], ncol=2 )
  return(merged)
}

#find clusters that should be merged based in covR (horizontal) distance of three clostest points
findMerge <- function(cluster, data, s, classes, min_distX, min_distY=0){
  
  merged <- matrix(NA, nrow=length(sel), ncol=2)
  
  #find smallest distance for each class pair
  for (r in seq_along(s)){
    clust1 <- data[cluster$cluster==classes$c1[s[r]] & ! is.na (cluster$cluster),]
    clust2 <- data[cluster$cluster==classes$c2[s[r]] & ! is.na (cluster$cluster),]
     
     # choose points of clusters that lie within 25% and 75% dh quantiles
     quan1 = quantile( clust1[, 2], probs=c(0.25,0.75) )
     quan2 = quantile( clust2[, 2], probs=c(0.25,0.75) )
     sub1 <- which( clust1[, 2] < quan1[2] & clust1[, 2] > quan1[1] )
     sub2 <- which( clust2[, 2] < quan2[2] & clust2[, 2] > quan2[1] )
   
    # get cov R distance and merge if smaller min_dist
    if (length(sub1) > 2 & length(sub2) > 2){
      dists   <- findMinimum( dist2(clust1[sub1,1], clust2[sub2,1], method='manhattan') )
      distsX  <- dists[,3]

      if ( sum(distsX < min_distX) == 3 ){
        merged[r,] <- c( max(classes[s[r], 1:2]), min(classes[s[r], 1:2]) ) 
      }
    }
  }
  
  merged <- matrix(merged[ ! is.na(merged) ], ncol=2 )
  return(merged)
}

#merge clusters defined in rows of merged and sum up memberships for merged classes
mergeCluster <- function(members, merged){
  
  mergeO <- matrix(merged[order(merged[,1], merged[,2], decreasing=TRUE ),], ncol=2)
  i <- 1
  
  while (nrow(mergeO) >0  ){
    maxClust <- max(mergeO[i,])
    minClust <- min(mergeO[i,])
    
    members[,minClust] = members[,minClust] + members[,maxClust]
    members <- members[, -maxClust]
    mergeO <- matrix(mergeO[-i,], ncol=2)
    
    #replace maxClust number in remaining matrix
    if (nrow(mergeO)> 0 & any(mergeO==maxClust)){
      ord <- which(rowSums(mergeO==maxClust) > 0)
      mergeO[mergeO==maxClust] = minClust
      
      #class with higher number should always be in first column of matrix 
      for (r in ord){
        mergeO[r,] <- sort(mergeO[r,],decreasing=TRUE)
      }
      
      #remove duplicates in table created by replacement
      mergeO <- matrix(mergeO[! duplicated(mergeO),], ncol=2)
      mergeO <- matrix(mergeO[order(mergeO[,1], mergeO[,2], decreasing=TRUE ),], ncol=2)
    }
    
  }
  return(as.matrix(members))
}

findMinimum <- function(dists){
  minima <- matrix(0, nrow=3, ncol=3)
  for (i in 1:3){
    minima[i,1:2] <- which(min(dists)==dists, arr.ind=TRUE)[1,]
    minima[i,3]   <- min(dists)
    dists[minima[i,1],] <- 100
    dists[,minima[i,2]] <- 100
  }
  return(minima)
}


removeOutlierEllipse <-function(cluster){

  coVar <- apply(cluster$parameters$variance$sigma, 3, function(i) i) 
  rem <- c()
  if (any(coVar[,,1][,1]>1)){
    rem <- which(coVar[,,1][,1]>1)
  }
  if (any(coVar[,,2][,2]>1)){
    rem <- c( rem, which(coVar[,,2][,2]>1) )
  }
  if (length(rem)>0)
    rem <- sort(unique(rem), decreasing=T)
  
  #before I used the cluster which was larger than 1 sd of dh values, which is not good for samples with spread dh values
  
  if (length(rem) >0){
    # remove all data for clusters that correspond to big ellipses
    cluster$parameters$mean <- cluster$parameters$mean[,-rem]
    cluster$z <- cluster$z[,-rem]
    ellipseX <- coVar[,,1][-rem,] #variance in X direction
    ellipseY <- coVar[,,2][-rem,] #variance in Y direction

    # set class to NA for those clusters that are removed
    for (out in rem){
      cluster$classification[cluster$classification==out] = NA
      sel <- cluster$classification > out & ! is.na(cluster$classification)
      cluster$classification[sel] = cluster$classification[ sel ]-1
    }

  }else{	
    ellipseX <- coVar[,,1] #variance in X direction
    ellipseY <- coVar[,,2] #variance in Y direction
  }
  return(list(cluster=cluster, ellipseX=ellipseX[,1], ellipseY=ellipseY[,2]))
}


removeOutlierPoints <- function(data, cluster_mat){
  cluster  <- data$cluster
  ellipseX <- data$ellipseX
  ellipseY <- data$ellipseY
  
  for (cl in unique(cluster$classification)){
    if (is.na(cl)){
      next
    }
    p <- which(cluster$classification==cl)
    #take median of standard deviation as allowed distance for closest points
    max_distX = 2*median(sqrt(ellipseX))
    max_distY = 2*median(sqrt(ellipseY)) #new
    
    dists <- apply(dist2(cluster_mat[p,],cluster_mat[p,],method='euclidean'), 1, function(i) sort(i)[2:4])
    #in case there are less than 3 points in the cluster
    dists <- data.frame(dists[! is.na(dists[,1]),])
    if (nrow(dists)>0){
    sel <- unlist(sapply(1:ncol(dists),function(i) if(any(dists[,i]>sqrt(max_distX^2+max_distY^2))) i))

    if(length(sel)>0)
      cluster$classification[p[unique(sel)]] <- NA
    }
  }
	
  return(cluster)
}

# 
removeOutlierPoints_cmean <- function(data, cluster_mat, massCenter, mainCluster){
  #get median range of points per cluster
  max_distX <- 0.7*median(sapply(sort(unique(data$cluster)), function(i){
                    s <- which(CM$cluster==i)
                    quantile(cluster_mat[s,1], p=0.95) - quantile(cluster_mat[s,1],p=0.05)
                }))
  max_distY <- 0.7*median(sapply(sort(unique(data$cluster)), function(i){
    s <- which(CM$cluster==i)
    quantile(cluster_mat[s,2], p=0.95) - quantile(cluster_mat[s,2],p=0.05)
  }))
  #max_distX = abs(quantile(cluster_mat[tmp,1], p=0.95) - quantile(cluster_mat[tmp,1],p=0.05))
  #max_distY = 1.2*abs(quantile(cluster_mat[tmp,2], p=0.95) - quantile(cluster_mat[tmp,2],p=0.05))
  
  for (cl in unique(data$cluster)){
    if (is.na(cl)){
      next
    }
    p <- which(data$cluster==cl)
    #take median of standard deviation as allowed distance for closest point

    x_dist <- abs(cluster_mat[p,1]-massCenter[cl,1])
    y_dist <- abs(cluster_mat[p,2]-massCenter[cl,2])
    
    sel <- which( x_dist > max_distX | y_dist > max_distY )
    #dont remove points from main cluster in X direction as these can't be estimated better than FWHM
    if( cl == mainCluster){
      sel <- which(y_dist > max_distY)
    }
    if(length(sel)>0)
       # points(cluster_mat[p[sel],], col='blue', pch=20)  
        data$cluster[p[sel]] <- NA
        if (length(sel) == length(p)){
          data$centers[cl,] <- NA
        }
  }
  
  return(data)
}


plotCov <- function(seg, chrLen){
  xtotal  = chrLen/ 10
  len     = chrLen/1000000
  # scale values
  seg$start       <- (seg$start/10)/xtotal
  seg$end         <- (seg$end/10)/xtotal
  
  sel <- which(seg$length < 1e6)
  
  # SNPs as data points, colored according to gain, loss and neutral
  
  p <- ggplot(environment=environment())  
  p <- p + geom_segment( data=seg, aes( x=start,y=tcnMean, xend=end, yend=tcnMean, col=as.character(cluster)), size = 1 )
  if (length(sel) >0)
    p <- p + geom_errorbar(data = seg[sel,], aes(x=start, ymin=tcnMean-0.01, ymax=tcnMean+0.01, col=as.character(cluster)))
  p <- p + scale_color_manual(values=c(col[1:length(unique(seg$cluster))], "grey"), name="cluster")
  # labs, title, boundaries
  p <- p + theme( legend.position="none", panel.grid=element_blank() )
  p <- p + xlab('Genomic position') + ylab('coverage') + xlim( c(0,1) ) +ylim(c(0,2))
  p <- p + theme( title = element_text(size=15), axis.title = element_text(size=12) )
  
  p1 <- ggplot(environment=environment())  
  p1 <- p1 + geom_segment( data=seg, aes( x=start,y=dhMax, xend=end, yend=dhMax, col=as.character(cluster)), size = 1 )
  if (length(sel)>0)
    p1 <- p1 + geom_errorbar(data = seg[sel,], aes(x=start, ymin=dhMax-0.01, ymax=dhMax+0.01, col=as.character(cluster)))
  p1 <- p1 + scale_color_manual(values=c(col[1:length(unique(seg$cluster))], "grey"), name="cluster")
  p1 <- p1 + xlab('Genomic position') + ylab('dhMax') + xlim( c(0,1) ) +ylim(c(0,1))
  p1 <- p1 + theme( legend.position="none", panel.grid=element_blank() )
  
  return(list(p, p1)) 
}


categorize <- function(data, min, max, intervall){
  categories <- data
  for( i in seq(max, min, -intervall) ){
    sel <- data < i
    categories[sel] <- paste0("<", i)
  }
  sel <- data>=max
  categories[sel] <- paste0(">", max)
  return(categories)
}


mergeClusters <- function(CM.tmp, minTcnMean, maxTcnMean, cluster_mat, mainCluster, minIncluded=0.85, MADfactor=3){
  # check how many points per cluster are within FWHM range around main Cluster, select those with >95% included
  identical = FALSE
  nCluster <- length(unique(CM.tmp$cluster))
  repeat{
    if( ! identical ){
  
      included <- sapply(sort(unique(CM.tmp$cluster)), function(cluster) {
        subCluster <- which(CM.tmp$cluster==cluster)
        sum( minTcnMean < cluster_mat[subCluster,1] & cluster_mat[subCluster,1] < maxTcnMean)/length(subCluster)
      })
      sel <- which(included >= minIncluded)
      if (length(sel) > 0 ){
        selName <- sort(unique(CM.tmp$cluster))[sel]
        #estimate median absolute deviation for these clusters and check whether any of them is broader than 
        med <- median(cluster_mat[CM.tmp$cluster==mainCluster,2])
        MADmain <- median( abs(cluster_mat[CM.tmp$cluster==mainCluster,2] - med)  )
        
        selMerge <- which( abs(CM.tmp$center[mainCluster,2] - CM.tmp$center[selName,2]) < MADfactor*MADmain )
        cat("merge ", names(selMerge), "\n")
        if(length(selMerge)>0){
          for (cl in selName[selMerge]){
            CM.tmp$cluster[CM.tmp$cluster==cl] <- mainCluster
            CM.tmp$centers[cl,] <- NA
          }
        }
      }
      
      CM.tmp$centers[mainCluster,1] <- mean(cluster_mat[CM.tmp$cluster==mainCluster,1]) 
      CM.tmp$centers[mainCluster,2] <- mean(cluster_mat[CM.tmp$cluster==mainCluster,2])
      
      if (length(unique(CM.tmp$cluster)) == nCluster){
        identical = TRUE
      }else{
        nCluster = length(unique(CM.tmp$cluster))
      }
    }else
      break
  }
  return(CM.tmp)
  
}

mergeClustersWithoutOutlier <- function(seg, leftCov, rightCov, center, mainCluster, minIncluded=0.85, MADfactor=3){
  # check how many points per cluster are within FWHM range around main Cluster, select those with >95% included
  identical = FALSE
  nCluster <- length(unique(seg$cluster[! is.na(seg$cluster)]))
  repeat{
    if( ! identical ){
  
      included <- sapply(sort(unique(seg$cluster[! is.na(seg$cluster)])), function(cluster) {
        subCluster <- which(seg$cluster==cluster)
        sum( leftCov < seg$tcnMean[subCluster] & seg$tcnMean[subCluster] < rightCov)/length(subCluster)
      })
      sel <- which(included >= minIncluded)
      if (length(sel) > 0 ){
        selName <- sort(unique(seg$cluster[! is.na(seg$cluster)]))[sel]
        #estimate median absolute deviation for these clusters and check whether any of them is broader than 
        med <- median(seg$dhMax[which(seg$cluster==mainCluster)])
        MADmain <- median( abs(seg$dhMax[which(seg$cluster==mainCluster)] - med)  )
      
        selMerge <- which( abs(center[mainCluster,2] - center[selName,2]) < MADfactor*MADmain )
        cat("merge ", rownames(center[selMerge,]), "\n")
        if(length(selMerge)>0){
          for (cl in selName[selMerge]){
            seg$cluster[which(seg$cluster==cl)] <- mainCluster
            center[cl,] <- NA
          }
        }
      }

      center[mainCluster,1] <- mean( seg$tcnMean[which(seg$cluster==mainCluster)] ) 
      center[mainCluster,2] <- mean( seg$dhMax[which(seg$cluster==mainCluster)] )
      
      if (length(unique(seg$cluster[! is.na(seg$cluster)])) == nCluster){
        identical = TRUE
      }else{
        nCluster = length(unique(seg$cluster[! is.na(seg$cluster)]))
      }
    }else
      break
  }
  return(list(seg,center))
  
}

mergePointsBySNP <- function(segments, clusterCenter, mainCluster, covValueRight, covValueLeft, minNumHets=5){
  #merge tose points in with less than minNumHets heterozygous SNPs in segment and that lie within 2 FWHM around center 
  sel <- which(segments$neighbour=="identical" & segments$tcnNbrOfHets < minNumHets &
    segments$tcnMean < (covValueRight) & segments$tcnMean > (covValueLeft) &
    is.na(segments$SV.Type) & ! is.na(segments$cluster))
  if(length(sel)>0)
    segments$cluster[sel] <- mainCluster
  return(segments)
  
}

mergePointsByError <- function( segments, clusterCenter, mainCluster ){
  #merge tose points in with less than minNumHets heterozygous SNPs in segment and that lie within 2 FWHM around center 
  sel <- which(segments$neighbour=="identical" & 
    ( segments$errorLength > segments$distTcn & segments$errorSNP > segments$distDH ) &
    is.na(segments$SV.Type) & ! is.na(segments$cluster))
  
  if (length(sel)>0)
    segments$cluster[sel] <- mainCluster
  return(segments)
  
}

findNeighbours <- function(segments, mainCluster){
  neighbour <- sapply(seq_len(nrow(segments)), function(i){
    if ( i==1 | i==nrow(segments) | any(segments$chromosome[(i-1):(i+1)] != segments$chromosome[i]) ){
      "edge"
    }else if ( is.na(segments$cluster[i]) | is.na(segments$cluster[i-1]) | is.na(segments$cluster[i+1]) ){
      NA
    }else if ( ! is.na(segments$cluster[i]) & segments$cluster[i]== mainCluster ){
      "self"
    }else if  ( ( ( segments$start[i] == segments$end[i-1] | segments$start[i] == segments$end[i-1]+1) & segments$cluster[i-1]== mainCluster ) &&
                  ( segments$end[i] == segments$start[i+1] |  segments$end[i] == segments$start[i+1]-1 ) & segments$cluster[i+1]== mainCluster ){
      "identical"
    }else
      NA
  } )
  return(neighbour)
}


combineNeighbours <- function(segments){
  segmentList <- split(segments, segments$chromosome)
  segmentList <- lapply(segmentList, function(chr)  {
    sel <- seq_len( nrow(chr) )
    for (j in seq_len(nrow(chr))) {
      repeat{
        if (  j < nrow(chr) && 
                chr$map[ sel[j] ] != 'homozygousDel' &&
                !is.na(chr$cluster[ sel[j] ]) &&
                chr$cluster[ sel[j] ] != "NA" &&
                !is.na(chr$cluster[ sel[j] + 1 ]) &&
                chr$cluster[ sel[j] ] == chr$cluster[ sel[j+1] ] &&
                is.na(chr$SV.Type[ sel[j]+1 ]) &&
                ( chr$end[ sel[j] ] == chr$start[sel[j+1]] | chr$end[sel[j]]+1 == chr$start[sel[j+1]] ) ) { 
          
          cat(paste0("",chr$start[sel[j]], " ",chr$end[sel[j]], " -> prune -> (",chr$start[sel[j]],":",chr$end[ sel[j + 1 ] ],")\n\n"))
          
          chr$end[ sel[j] ] = chr$end[ sel[j + 1 ] ]
          chr$length[ sel[j] ] = chr$length[ sel[j] ] + chr$length[ sel[j + 1] ]
          chr$tcnMean[ sel[j] ] = (chr$tcnNbrOfLoci[ sel[j] ]*chr$tcnMean[ sel[j] ] + chr$tcnNbrOfLoci[ sel[j + 1] ]*chr$tcnMean[ sel[j + 1] ]) / 
            (chr$tcnNbrOfLoci[ sel[j]] + chr$tcnNbrOfLoci[ sel[j + 1] ])
          
          chr$dhMax[ sel[j] ] = (chr$tcnNbrOfHets[ sel[j] ]*chr$dhMax[ sel[j] ] + chr$tcnNbrOfHets[ sel[j + 1] ]*chr$dhMax[ sel[j + 1] ]) /
            (chr$tcnNbrOfHets[ sel[j] ] + chr$tcnNbrOfHets[ sel[j + 1] ] )
          
          chr$tcnNbrOfLoci[ sel[j] ] = chr$tcnNbrOfLoci[ sel[j] ] + chr$tcnNbrOfLoci[ sel[j + 1] ]
          chr$tcnNbrOfSNPs[ sel[j] ] = chr$tcnNbrOfLoci[sel[j]]
          
          chr$dhNbrOfLoci[ sel[j] ] = chr$dhNbrOfLoci[ sel[j] ] + chr$dhNbrOfLoci[ sel[j + 1] ] #This is not used at a later stage. Therefore it is not important to change this (original calculation not known)
          chr$tcnNbrOfHets[ sel[j] ] = chr$tcnNbrOfHets[ sel[j] ] + chr$tcnNbrOfHets[ sel[j + 1] ] #This is used later and thus adjusted properly to the definition
          chr = chr[-(sel[j + 1]), ]
          sel[ sel > sel[j + 1] ] = sel[ sel > sel[j + 1] ]-1
          sel = sel[-(j + 1)]
        } else {
          break
        }
      }
    }
    chr
  })
  return( do.call(rbind, segmentList))
}  


generatePlots<- function(segments, selIdentical, mainCenter, covLeftHalf, covRightHalf, covLeftFull, covRightFull, minNbrOfHets = 5){

     pSnpNbr <- ggplot(data=segments, aes(log2(tcnMean), dhMax) ) +
            geom_point() + 
       	    ggtitle("nbrOfHets< minNbrOfHets" ) + 
            scale_color_manual(values=c('red','blue'), name='nbrOfHets<5') +
            geom_vline( xintercept=c( log2(covLeftFull) ), col='red')+
            geom_vline( xintercept=c( log2(covRightFull) ), col='red' ) +
            geom_vline( xintercept=c( log2(covLeftHalf) ), col='blue')+
            geom_vline( xintercept=c( log2(covRightHalf) ), col='blue' ) +
            theme(legend.position="none")

  pSNPError <- ggplot(data=segments, aes(log2(tcnMean), dhMax) ) +
            geom_point() + 
	    ggtitle("errorSNP > distDH") + 
            scale_color_manual(values=c('red','blue'), name='errorSNP>distDH') +
            geom_vline( xintercept=c( log2(covLeftFull) ), col='red')+
            geom_vline( xintercept=c( log2(covRightFull) ), col='red' ) +
            geom_vline( xintercept=c( log2(covLeftHalf) ), col='blue')+
            geom_vline( xintercept=c( log2(covRightHalf) ), col='blue' ) +
            theme(legend.position="none")

  pError <- ggplot(data=segments, aes(log2(tcnMean), dhMax) ) +
            geom_point() + 
	    ggtitle("errorLength >distTcn & errorSNP > distDH")  + 
            scale_color_manual(values=c('red','blue'), name='totalErr>totalDist' ) +
            geom_vline( xintercept=c( log2(covLeftFull) ), col='red')+
            geom_vline( xintercept=c( log2(covRightFull) ), col='red' ) +
            geom_vline( xintercept=c( log2(covLeftHalf) ), col='blue')+
            geom_vline( xintercept=c( log2(covRightHalf) ), col='blue' ) +
            theme(legend.position="none")


  pLength <- ggplot(data=segments, aes(log2(tcnMean), dhMax) ) +
            geom_point() + 
	    ggtitle( "errorLength>distTcn" ) + 
            scale_color_manual(values=c('red','blue'), name='neighbourMain' ) +
            geom_vline( xintercept=c( log2(covLeftFull) ), col='red')+
            geom_vline( xintercept=c( log2(covRightFull) ), col='red' ) +
            geom_vline( xintercept=c( log2(covLeftHalf) ), col='blue')+
            geom_vline( xintercept=c( log2(covRightHalf) ), col='blue' ) +
            theme(legend.position="none")

  if ( length(selIdentical) > 0 ){
	pSnpNbr   <- pSnpNbr + geom_point(data=segments[selIdentical,], aes(x=log2(tcnMean), y=dhMax, col=( tcnNbrOfHets < minNbrOfHets) ) )
	pSNPError <- pSNPError + geom_point(data=segments[selIdentical,], aes(x=log2(tcnMean), y=dhMax, col=( errorSNP >distDH) ) )
	pError	  <- pError  + geom_point( data=segments[selIdentical,], aes(x=log2(tcnMean), y=dhMax, col=( (errorLength > distTcn & errorSNP > distDH)  ) ) ) 
	pLength	  <- pLength + geom_point(data=segments[selIdentical,], aes(x=log2(tcnMean), y=dhMax, col=(errorLength>distTcn) ) ) 
  }


  p <- grid.arrange(pSnpNbr, pSNPError, pError, pLength,ncol=2)
  return(p)
}


mergePoints <- function(segments, mainCenter, mainCluster, covLeftValue, covRightValue, minNumHets=5){
  i=0
  segments$distDH <- abs(segments$dhMax -mainCenter$dhMax)
  segments$errorSNP <- 1/sqrt(segments$tcnNbrOfHets)
  segments$distTcn <-  abs( segments$tcnMean - mainCenter$tcnMean )
  segments$errorLength  <- 2/log2(segments$length)
  segments$totalError <- segments$errorLength +segments$errorSNP
  identical <- FALSE
  cluster_old <- length(segments$cluster)
  
  repeat{
    if( ! identical ){
      segments <- mergePointsBySNP(segments, mainCenter, mainCluster, covLeftValue, covRightValue, minNumHets=minNumHets)
      segments <- combineNeighbours(segments)
      segments <- mergePointsByError(segments, mainCenter, maxCluster)
      segments$neighbour <- findNeighbours(segments, mainCluster)
      segments <- combineNeighbours(segments)
      segments$neighbour <- findNeighbours(segments, mainCluster)
      
      segments$distDH <- abs(segments$dhMax -mainCenter$dhMax)
      segments$errorSNP <- 1/segments$tcnNbrOfHets
      segments$distTcn <-  abs(segments$tcnMean -mainCenter$tcnMean)
      segments$errorLength  <- 2/log2(segments$length)
      segments$totalError <- segments$errorLength +segments$errorSNP
      if( cluster_old == length(segments$cluster) )
        identical=TRUE
      cluster_old <- length(segments$cluster)
      i=i+1
    }else
      break
  }
  cat(i-1, " iterations completed\n")
  return(segments)
}


removeOutlierPoints_cmean_alt <- function(data, center, deviationFactor = 2){
  
  #determine distance of each point to corresponding center
  rowMin <- replicate(nrow(data), NA)
  
  for (cl in as.numeric( unique( data$cluster[ data$cluster!="NA" ]) ) ){
    if ( is.na(cl) | cl == "NA"){
      next
    }
    p <- which(data$cluster==cl)
    #take median of standard deviation as allowed distance for closest point
    allDist <- dist2( data[p, c("tcnMean", "dhMax")], center[as.character(cl), ] )
    rowMin[p]  <- allDist[,1]
  }
  
  rowMinAverage <- rowMin/mean(rowMin,na.rm=T)
  #remove those that deviate more that determined factor
  sel <- rowMinAverage > deviationFactor  | is.na(rowMinAverage)
  #main <- which(data$cluster == mainCluster)            
  #sel[main] <- FALSE
  
  ggplot( data.frame(data), aes(log2(tcnMean), dhMax, col=sel)) +#, col=sel ) )  + 
    geom_point() + geom_point(data=data.frame(center), aes(log2(tcnMean), dhMax), col='red' )
  
  data$cluster[sel] <- NA
  
  return(data)
}
