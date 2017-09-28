

# Copyright (c) 2017 The ACEseq workflow developers.
# This script is licenced under (license terms are at
# https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
################################################################################
## functions for purity_ploidy.R
################################################################################
library(flux)

# == title
# computing by chromosomes
#
# == param
# -chromosomes chromosome name
#
# == value
# character vector
runTheStuff = function(segments, chromosomes, minLim, maxLim, averageCov, minCov) {

  index = 1
	k = lapply(chromosomes, function(chr) doItForOneChr(segments,chr, minLim, maxLim, averageCov, minCov))

	myPeaks2=c()
	for(i in seq_along(k)) {
		for (j in seq_along(k[[i]])) {
			myPeaks2 = c(myPeaks2, k[[i]][[j]])
		}
	}
	
	return(myPeaks2)
}



# == title
# calculates area under curve left and right from 0.5
#
# == param
# -chr chromosome index
#
# == value
# a vector
doItForOneChr = function(segments, chr, minLim, maxLim, averageCov, minCov) {
  cat("chromosome: ", chr, "\n")
	## variables from parent environment:
	## segments
	## dataAll

	allAreas = c()
	
	sel = which(segments$chromosome == chr)
	start = segments$start[sel]

	index = 1
	
	## save images in a multiple page PDF
	pdf(paste0(out,"/chr_", chr,"_peaks.pdf"), width = 5, height = 5)
	
	for (seg in start) {
	
		## add 2013-9-25
		cat(paste0( "segment: ", seg, "\n\n"))
		
		sel <- which(dataAll[[chr]]$start == seg &
		             dataAll[[chr]]$betaN > 0.3 & 
		             dataAll[[chr]]$betaN < 0.7 & 
		             ! is.na(dataAll[[chr]]$betaT) & 
		             #dataAll[[chr]]$Atumor  + dataAll[[chr]]$Btumor > 20  & # removed 
		             #dataAll[[chr]]$Atumor  + dataAll[[chr]]$Btumor < 80  & # removed
		             dataAll[[chr]]$Anormal + dataAll[[chr]]$Bnormal > minCov & 
		             dataAll[[chr]]$Anormal + dataAll[[chr]]$Bnormal < averageCov*2.5 ) 
		
      if (length(sel) > 30) {
      #draw density
	tmp = density(dataAll[[chr]]$betaT[sel], bw=0.05)
			#see ratio of area
	limit <- 0.5
	temp_start <- min(tmp$x)
	temp_stop <- max(tmp$x)
	left_ind <- which(tmp$x <= limit)
  if (length(left_ind)>1){
	  area_left <- auc(tmp$x[left_ind],tmp$y[left_ind])
  }else{
    area_left <- 0
  }
	right_ind <- which(tmp$x > limit)
  if (length(right_ind) > 1){
	  area_right <- auc(tmp$x[right_ind],tmp$y[right_ind])
  }else{
    area_right <- 0
  }
	position_quot <- area_right / area_left+0.000001
	position_rel_diff <- abs(area_right - area_left) / (area_right + area_left)
	# can also potentially include an asymmetry criterion here
	plot(tmp, col="blue", main = paste0("chr: ",chr, "  segment: ", index,", SNPS: ",length(sel),"\narea_left = ",
						round(area_left,digits=3),"; area_right = ",round(area_right, digits=3),
						"\nposition_quot = ", round(position_quot,digits=3),"; position_rel_diff = ",
						round(position_rel_diff, digits=3)), xlab="BAF", xlim = c(-0.2, 1.2),
	      sub=paste0("start: i",seg))               
	polygon(c(tmp$x[left_ind],limit),c(tmp$y[left_ind],0),col="red")
	polygon(c(limit,tmp$x[right_ind]),c(0,tmp$y[right_ind]),col="blue")
	abline(v=0.5) 
      
	
	allAreas = c(allAreas, position_rel_diff)
      }else{
	allAreas = c(allAreas, NA)
      }
      index = index + 1
      }
	
      dev.off() ##
      cat(chr, "\t", length(allAreas), "\t", length(start), "\n")
      return(allAreas)
}
