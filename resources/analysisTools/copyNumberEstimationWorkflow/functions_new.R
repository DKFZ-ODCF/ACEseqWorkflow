
################################################################################
## functions for purity_ploidy.R
################################################################################

library(multicore)

# == title
# calculate the derivative of a curve, I guess
# 
# ==param
# -x x coordinate
# -y y coordinate
# -width
# -sigma
#
# ==value
# a data frame with coordinate of derivative curve
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

# == title
# find peaks for one chromosome
#
# == param
# -chr chromosome index
#
# == value
# a vector
doItForOneChr = function(chr) {

	## variables from parent environment:
	## segments
	## dataAll

	allPeaks = c()
	
	sel = which(segments$chromosome == chr)
	start = segments$start[sel]

	index = 1
	
	## save images in a multiple page PDF
	pdf(qq("@{out}/chr_@{chr}_peaks.pdf"), width = 5, height = 5)
	
	for (seg in start) {
	
		## add 2013-9-25
		cat(qq("segment: @{seg}\n\n"))
		
		sel <- which(dataAll$chromosome == chr & 
		             dataAll$start == seg &
		             dataAll$betaN > 0.3 & 
		             dataAll$betaN < 0.7 & 
		             dataAll$Atumor + dataAll$Btumor > 20 & 
		             dataAll$Atumor + dataAll$Btumor < 80 & 
		             dataAll$Anormal + dataAll$Bnormal > 20 & 
		             dataAll$Anormal + dataAll$Bnormal < 80)
		
		if (length(sel) > 30) {
			tmp = density(dataAll$betaT[sel], bw=0.05)
			dev = derive(tmp$x,tmp$y,width=25,sigma=9)  #width was 15 before

			zeroCrossings <- which(diff(sign(dev$y))!=0) # finds points where sign changes
			BAFs <- dev$x[zeroCrossings] # BAF value at the position where the sign changes
			center=which(BAFs > 0.45 & BAFs < 0.55)

			minimaInCross = zeroCrossings[zeroCrossings %in% which(diff(sign(dev$y)) > 0)]# positive sign -> minima
			maximaInCross = zeroCrossings[zeroCrossings %in% which(diff(sign(dev$y)) < 0)] # negative sign -> maxima

			minima = cbind(dev$x[minimaInCross], tmp$y[minimaInCross])
			minima = data.frame(minima)
			colnames(minima) = c("BAF","density")

			maxima = cbind(dev$x[maximaInCross], tmp$y[maximaInCross])
			maxima = data.frame(maxima)
			colnames(maxima) = c("BAF", "density")

			valid_minima = minima[minima$BAF > 0.45 & minima$BAF < 0.55, , drop = FALSE]  ## add drop, 2013-9-25
			valid_maxima = maxima[maxima$BAF > 0.45 & maxima$BAF < 0.55, , drop = FALSE]  ## add drop, 2013-9-25
			
			peaks = NA
			reason = NA
			if (nrow(valid_minima) > 0) { #has the density a minimum close to 0.5?
			
				if (nrow(valid_minima) == 1) { #only one minimum between 0.45 & 0.55?
					#identify maxima left and right
					maxima_left = maxima[maxima$BAF < valid_minima$BAF, , drop = FALSE]   ## add drop, 2013-9-25
					maxima_right = maxima[maxima$BAF > valid_minima$BAF, , drop = FALSE]   ## add drop, 2013-9-25
	
					## original code is abstracted as a function
					res = evaluate_maxima(maxima_left, maxima_right)
					peaks = res$peaks
					reason = res$reason
					
				} else { #more than one minimum between 0.45 & 0.55?
					
					peaks = "crap"
					reason = "more than one maximum between 0.45 and 0.55"
				}
			} else { # no minimum close to 0.5
				# hier sollte wenn dann nur noch ein maximum im valid bereich liegen

				if (nrow(valid_maxima) == 1) { #gibt es ein maximum zwischen 0.45 und 0.55
					if (valid_maxima$density >= max(maxima$density)) {
					    #maximum auch globales maximum
						peaks = 1
						reason = "one central maximum"
					} else { #maximum nicht globales maximum
						maxima_left = maxima[maxima$BAF < valid_maxima$BAF, , drop = FALSE]   ## add drop, 2013-9-25
						maxima_right = maxima[maxima$BAF > valid_maxima$BAF, , drop = FALSE]   ## add drop, 2013-9-25
						
						## original code is abstracted as a function
						res = evaluate_maxima(maxima_left, maxima_right)
						peaks = res$peaks
						reason = res$reason
					}

				} else { #kein maximum im valid bereich
					maxima_left = maxima[maxima$BAF < 0.5, , drop = FALSE]   ## add drop, 2013-9-25
					maxima_right = maxima[maxima$BAF > 0.5, , drop = FALSE]   ## add drop, 2013-9-25
	
					## original code is abstracted as a function
					res = evaluate_maxima(maxima_left, maxima_right)
					peaks = res$peaks
					reason = res$reason
				}
			}

		} else {
			peaks = "missing"
		}

		if (peaks == 2) {
		
			cat(qq("maxima_left: @{maxima_left}\n\n"))
			cat(qq("maxima_right: @{maxima_right}\n\n"))
			
			central_maxima_left = maxima_left[maxima_left$BAF > 0.45, , drop = FALSE]
			central_maxima_right = maxima_right[maxima_right$BAF < 0.55, , drop = FALSE]
			
			if (length(central_maxima_left) > 0 && 
			    length(central_maxima_right) > 0 && 
			    (max(maxima$density) %in% central_maxima_left$density || max(maxima$density) %in% central_maxima_right$density)) {
				
				peaks=1
				reason="combined central maximum"
			}
		}
		
		# plots into multiple pages PDF
		if (peaks != "missing") {   
			
			plot(tmp, col="blue", main = qq("chr: @{chr}, segment: @{index}, peaks: @{peaks},\nreason: @{reason} SNPS:@{length(sel)}"), xlab="BAF", xlim = c(-0.2, 1.2))               
			lines(dev$x, (dev$y * 10) + 1, col = "red")
			lines(dev$x, (dev$y * 10) * 0 + 1, col = "red", lty = 2)
			abline( v=c(0.45,0.55), lty= 2 )
			abline( v=0.5, lty=3 )
			
		}

		allPeaks = c(allPeaks, peaks)
		index = index + 1
	}
	
	dev.off() ##
	#print(allPeaks)
	return(allPeaks)
}


# == title
# parallel computing by chromosomes
#
# == param
# -chromosomes chromosome name
# -threads     number of threads
#
# == value
# character vector
runTheStuff = function(chromosomes = 1:5, threads = 12) {

	index = 1
	k = list()
	jobs = list()
	jobsIndex = 0
	
	for (chr in chromosomes) {
	
		jobsIndex = jobsIndex+1
		jobs[[jobsIndex]] = parallel(doItForOneChr(chr), silent = FALSE)
		if(jobsIndex == threads || chr == chromosomes[length(chromosomes)]) {
			k[[index]] = collect(jobs, wait = TRUE)
			index = index + 1
			jobsIndex = 0
		}
	}

	myPeaks2=c()
	for(i in seq_along(k)) {
		for (j in seq_along(k[[i]])) {
			myPeaks2 = c(myPeaks2, k[[i]][[j]])
		}
	}
	
	return(myPeaks2)
}

## extract from `doItForOneChr`
evaluate_maxima = function(maxima_left, maxima_right) {


	heights = c(maxima_left$density, maxima_right$density)
	heights_rat = sort( heights/max(heights), decreasing = TRUE )
  if ( length(heights_rat) == 1 ) {
    peaks = 2
    reason = "Single maximum outside 0.45 and 0.55"
	}else if ( length(heights_rat) > 1 & heights_rat[2] < 0.05  ){
		peaks = 2
		reason = "Maximum outside 0.45 and 0.55"
	} else  if (nrow(maxima_left) > 0 && nrow(maxima_right) > 0) { #are there maxima left and right of the minimum
		 
		takeIt = c()
		for (i in seq_len(length(maxima_left) / 2)) {
			for (j in seq_len(length(maxima_right) / 2)) {
				if (abs(0.5 - maxima_left$BAF[i] - maxima_right$BAF[j] + 0.5) < 0.08) { #are the maxima symmetric?
					takeIt = rbind(takeIt, c(i, j))
				}
			}
		}
		
		if (length(takeIt) > 0) { #symmetric maxima found
				peaks = 2
				reason = "two symmetric maxima"
		} else{
				peaks = "crap"
				reason = "no symmetric maxima found"
		}

	}else{ #there are no maxima left and right of the minimum
		peaks = "crap"
		reason = "no maxima between 0.45 and 0.55 AND no symmetric maxima left and right"
	}
	
	return(list(peaks = peaks, reason = reason))
}
