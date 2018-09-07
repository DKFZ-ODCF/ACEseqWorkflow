##################

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
## start Daniel ##
##################

extractFWHM <- function(ratio_density){
  half_max_ratio_density <- max(ratio_density$y)/2
  max_pos <- ratio_density$x[which(ratio_density$y==max(ratio_density$y))]
  right_ind <- which(ratio_density$x>=max_pos)
  temp_diff_ratio_density <- abs(ratio_density$y-half_max_ratio_density)
  half_max_pos_right_ind <- which(temp_diff_ratio_density==min(temp_diff_ratio_density[right_ind]))
  half_max_pos_left_ind <- which(temp_diff_ratio_density==min(temp_diff_ratio_density[-right_ind]))
  half_max_pos_right <- ratio_density$x[half_max_pos_right_ind]
  half_max_pos_left <- ratio_density$x[half_max_pos_left_ind]
  FWHM <- half_max_pos_right-half_max_pos_left
  return(list(half_max_pos_right=half_max_pos_right,half_max_pos_left=half_max_pos_left,FWHM=FWHM))
}

################
## end Daniel ##
################

# rewrite the function defineMainCluster to also return the width of the main cluster for quantification and QC
defineMainCluster <- function(cov, fit, type="unknown", plotFlag=NULL){
  
  #Determine distribution of difference between fit and points
  my_diff			<- cov - fit$y
  my_diff_rel		<- my_diff/fit$y
  my_density		<- density(my_diff_rel, adjust=0.1,from=-1,to=10)
  
  if( !is.null(plotFlag) ){
    png(paste0(plotFlag,"/", type, "density_defineMainCluster.png"), type="cairo" )
      plot(my_density)
    dev.off()
  }
  ##comment: very noisy sample have several peaks at the maximum and lead to a very narrow main cluster
  ##due to a close minima, this won't be corrected as it will affect other samples
  ## e.g magic MDT-AP-2719, PC58s-000006
  
  #find absolute maximum of function and closest minima as borders of main cluster
  max_index = which(my_density$y==max(my_density$y))
  abs_max			<- my_density$x[ max_index ]
  my_density  <- remove_duplicated_values(my_density, max_index)
  my_local_minima_ind  	<- which(diff(sign(diff(my_density$y)))==2) + 1	#shift index by one due to double diff that removes first entry
  my_local_minima		    <- my_density$x[which(diff(sign(diff(my_density$y)))==2)]
  my_local_dist_from_max	<- my_local_minima - matrix(abs_max,length(my_local_minima))
  my_left_min_index	<- max(which(my_local_dist_from_max<0))
  my_left_min  	<- my_local_minima[my_left_min_index]
  if (my_left_min_index == -Inf){
    my_left_min_index <- 1
    my_left_min    <- my_density$x[1]
  }
  my_right_min_index	<- min(which(my_local_dist_from_max>0))
  my_right_min  	<- my_local_minima[my_right_min_index]
  if (my_right_min_index == Inf){
    my_right_min_index <- length(my_density(x))
    my_right_min <- my_density$x[length(my_density$x)]
  }
  
  #check the height of the minima compared to the maxima in order find extremly noisy samples
  abs_max_y <- max(my_density$y)
  diff_left <- abs_max_y - my_density$y[my_local_minima_ind[my_left_min_index] ]
  diff_right <- abs_max_y - my_density$y[my_local_minima_ind[my_right_min_index]]
  
  if(diff_left < 0.1*abs_max_y  | diff_right < 0.1*abs_max_y ){
    if (diff_left < 0.05*abs_max_y  | diff_right < 0.05*abs_max_y ){
      cat( "WARNING: One or both borders of the main cluster are poorly defined in ", type," ( red_flag ) \n") 
    }else{
      cat( "WARNING: One or both borders of the main cluster might be poorly defined in", type, " ( yellow_flag ) \n") 
    }
  }
  
  my_main_cluster_ind	<- which(my_diff_rel < my_right_min*scale_factor & my_diff_rel > my_left_min*scale_factor)

#   ##################
#   ## start Daniel ##
#   ##################
# 
#   #now compute FWHM of the cluster by calling other method
#   FWHM_data <- extractFWHM(my_density)  
#   return(list(ind=my_main_cluster_ind, width=my_right_min-my_left_min, FWHM=FWHM_data$FWHM, half_max_pos_right=FWHM_data$half_max_pos_right, half_max_pos_left=FWHM_data$half_max_pos_left))
#   
#   ################
#   ## end Daniel ##
#   ################

  return(list(ind=my_main_cluster_ind, width=my_right_min-my_left_min, dens=my_density))
}


remove_duplicated_values <- function(dens, max_index){
  dens_left  <- list( x=dens$x[1:max_index], y=dens$y[1:max_index] )
  dens_right <- list( x=dens$x[(max_index+1):length(dens$y)], y=dens$y[(max_index+1):length(dens$y)] )
  
  if ( any(diff(dens_left$y) == 0) ){
    #remove value further away from peak
    rem_left <- which(diff(dens_left$y)==0)
    dens_left$x <- dens_left$x[-rem_left]
    dens_left$y <- dens_left$y[-rem_left]
  }
  
  if ( any(diff(dens_right$y) == 0 ) ){
    #remove value furhter away from peak
    rem_right <- which(diff(dens_right$y)==0)+1
    dens_right$x <- dens_right$x[-rem_right]
    dens_right$y <- dens_right$y[-rem_right]
  }
  densNew <- list(x=c(dens_left$x, dens_right$x), y=c(dens_left$y, dens_right$y) )
                  
  return(densNew)
}

#function to get derivative of density curve by looking at "width" neighbouring points
derive <- function(x, y, width = 7, sigma = 1) {
  numOfWeights = width - 1
  weights = c()
  for (i in seq_len(numOfWeights / 2)) {
  	weights = c(weights, 1/(sigma*sqrt(2*pi))*exp(-0.5*(((i/sigma)^2))))
	}
	
	weights = c(rev(weights), weights)
	weights = weights * 1 / sum(weights)
	
	deriv = rep(NA, length(x))
	for (thisX in (numOfWeights/2+1):(length(x)-numOfWeights/2-1)) {
		deriv[thisX] = 0
		for (i in 1:numOfWeights) {
			deriv[thisX] = deriv[thisX] + weights[i] * (y[thisX + i - numOfWeights/2] - y[thisX + i - numOfWeights/2 - 1])
		}
	}
	derivative = cbind(x, deriv)
	
	colnames(derivative) = c("x", "y")
	return(as.data.frame(derivative))
}

#function to check whethetr density distribution over normalized coverage has one or two peaks and if
#the peak is further 
checkControl <- function(coverage, covIndex){

  par(mfrow=c(3,2))
  diffPeaks <- NULL
  for (chr in sort( unique( as.numeric(coverage$chromosome) ) ) ){
    cat(chr,"\n")
    if(chr>=23)
      next
    
    sel <- which(coverage$chromosome==chr)
    dens <- density(coverage[sel,covIndex])
    
    #get derivative
    dev <- derive(dens$x, dens$y, width=5,sigma=1)
    #find maxima    
    zeroCrossings <- which( diff( sign(dev$y) ) != 0) # finds points where sign changes
    maximaInCross = zeroCrossings[zeroCrossings %in% which(diff(sign(dev$y)) < 0)] 
    
    #check for second Peak
    maxPeak <- which(dens$y==max(dens$y[zeroCrossings])) 
    secondPeak <- maximaInCross[ which( dens$y[maximaInCross] >= 0.1*dens$y[maxPeak] & dens$y[maximaInCross] != dens$y[maxPeak] )[1] ]
    
    if (  0.5*( round(2*dens$x[maxPeak])) != 1  | ( ! is.na( secondPeak ) ) ) {
      cat( paste(chr, "warning indicator for contaminated sample or sample swap!\n") )
      if( is.na(secondPeak) ){
	      diffPeaks <- c(diffPeaks, "shifted")
      }else
        diffPeaks <- c( diffPeaks, abs(dens$x[maxPeak] - dens$x[secondPeak]) )

    }else{
        diffPeaks <- c(diffPeaks, NA)
    }

    plot(dens, col="blue", main = paste0("chr: ",chr), xlab="normalized coverage control")
  }
  
  return(diffPeaks)
}

#create coverage plots
plotCoverage <- function(coverageTab, chromosomeBorders=NULL, chr=NULL, ylims=c(-4,4) ){
    
    labelposition=NULL
    if ( ! is.null(chromosomeBorders) )
	labelPosition <- (diff(chromosomeBorders)/2+chromosomeBorders[1:24])/1e6
	
    par(mfrow=c(3,1), cex.lab=1.5, cex.main=2)
    
    plot(coverageTab$start/1000000, log2(coverageTab$covNnorm), xlab='genomic coordinate (MB)', ylab=' log2 corrected normalized coverage control', cex=0.01, main=paste0( "chr", chr, "\n", "Control coverage"), ylim=ylims )
    abline(h=0, col='red')
    if ( ! is.null(chromosomeBorders) ){
      abline( v=chromosomeBorders/1000000, lty='dotted' )
      text(x= labelPosition, y= 3.5, labels=as.character(c(1:22, 'X', 'Y')), cex=1.5 )
    }
    
    plot(coverageTab$start/1000000, log2(coverageTab$covTnorm), xlab='genomic coordinate (MB)', ylab='log2 corrected normalized coverage tumor', cex=0.01, main='Tumor coverage',  ylim=ylims )
    abline(h=0, col='red')
    if ( ! is.null(chromosomeBorders) ){
      abline( v=chromosomeBorders/1000000, lty='dotted' )
      text(x= labelPosition, y= 3.5, labels=as.character(c(1:22, 'X', 'Y')), cex=1.5 )
    }
    
    plot(coverageTab$start/1000000, log2(coverageTab$covR), xlab='genomic coordinate (MB)', ylab='log2 corrected normalized coverage ratio', cex=0.01, main='Tumor/Control coverage ratio',  ylim=ylims)
    abline(h=0, col='red')
    if ( ! is.null(chromosomeBorders) ){
      abline( v=chromosomeBorders/1000000, lty='dotted' )
      text(x= labelPosition, y= 3.5, labels=as.character(c(1:22, 'X', 'Y')), cex=1.5 )
    }
    
}

#function to adjust startcoordinates of chromosomes 
adjustCoordinates <- function(sortedLengthTab, coverageTab){

    addToStart <- 0
    chromosomeBorders <- 0
    
    for (chr in sortedLengthTab$chromosome ){
      sel <- which( coverageTab$chromosome== chr )
      if( length(sel) > 0 )
	coverageTab$start[sel] <- coverageTab$start[sel] + addToStart
	
      addToStart <- addToStart + sortedLengthTab$length[chr] 
      chromosomeBorders <- c(chromosomeBorders, addToStart)
    }  
    
  return( list(coverageTab=coverageTab, chromosomeBorders=chromosomeBorders) )
}
