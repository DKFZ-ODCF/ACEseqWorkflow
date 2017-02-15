require(ggbio)
require(ggplot2)
require(gridExtra)

#plot TCNs
#chr is the chromosome number
#ratio: data.frame; subset (chromosome) of dataAll dataframe containing all SNPs
#seg: data.frame; subset of combi which contains all segments for sample
plotTCN = function (chromosome, ratio, seg, Ploidy, Purity, fullPloidy, chrLen, ymaxcov, plots='single', crestSub=NULL, p = NULL) {
		
		ymaxDH	= 1 
		xtotal	= chrLen/ 10
		len	= chrLen/1000000

		# scale values
		seg$start    	<- (seg$start/10)/xtotal
		seg$end		<- (seg$end/10)/xtotal

		ratio$start	<- (ratio$start/10)/xtotal
		ratio$end	<- (ratio$end/10)/xtotal
		ratio$SNP	<- (ratio$SNP/10)/xtotal

		copyTColors = c("#000000","#228B22", "#8B0000")
		colScale <- scale_colour_manual( values = c( 'n'=copyTColors[1], 'g'=copyTColors[2], 'l'=copyTColors[3] ) )

		#limit plots to TCN 8 to avoid displaying high level amplifications
		if (ymaxcov>10){
			ymaxcov <- 8
		}

		# SNPs as data points, colored according to gain, loss and neutral
		p <- ggplot(environment=environment())
		if(nrow(ratio)>1){
		  p <- p + geom_point( data=ratio[ seq(1 ,nrow(ratio), 7), ], aes( x=SNP, y=copyT, col=GNL ), pch=16, cex=0.3 ) +
					colScale + theme_bw()
		}
		# allele specific and total copy numbers: c1Mean, c2Mean and tcnMean
	        p <- p + geom_segment( data=seg, aes( x=start, xend=end, y=c1Mean, yend=c1Mean ), colour="#00BFFFFF", size = 1, na.rm=TRUE )
		p <- p + geom_segment( data=seg, aes( x=start,xend=end, y=c2Mean, yend=c2Mean ), colour="#00BFFFFF", size = 1, na.rm =TRUE )
		p <- p + geom_segment( data=seg, aes( x=start,y=tcnMean, xend=end, yend=tcnMean ), colour="#0000CDFF", size = 1 )

		# labs, title, boundaries
		p <- p + theme( legend.position="none", panel.grid=element_blank() )
		p <- p + xlab('') + ylab('TCN') + xlim( c(0,1) )
    		p <- p + scale_y_continuous(limits=c(0, ymaxcov + 1.2), breaks=seq(0, ymaxcov+1.2, 2), labels=as.character( seq( 0, (ymaxcov+1.2), 2 ) ) )

		p <- p + theme( title = element_text(size=15), axis.title = element_text(size=12) ) 


		p <- p + geom_hline( yintercept=ymaxcov, col = "#000000FF", lty=1, lwd =0.2 )
		p <- p + geom_hline( yintercept=seq( 0, ymaxcov, 1 ),  col = "#C0C0C0", lty="dotted", lwd =0.5 )
		# allele specific and total copy numbers: c1Mean, c2Mean and tcnMean


		if ( any(seg$tcnMean>ymaxcov) ){
			highAmp <- which( seg$tcnMean > ymaxcov )
			highSeg <- seg[highAmp, ]
			highSeg$tcnMean <- 8
			highSeg$middle <- highSeg$start + ( highSeg$end - highSeg$start )/2
			p <- p + geom_segment( data=highSeg, aes( x=start, y=tcnMean, xend=end, yend=tcnMean ), colour="#0000CDFF", size = 1 )
			p <- p + geom_point( data=highSeg, aes( x=middle, y=tcnMean), colour="#0000CDFF", shape=17)
		}

		if (plots == 'single'){
			# segment boundaries
			p <- p + geom_vline( xintercept=seg$start, col = "#000000DD", lty=5, lwd =0.2 )
      			p <- p + geom_vline( xintercept=seg$end, col = "#000000DD", lty=5, lwd =0.2 )
			# add axis
			p <- p + scale_x_continuous( breaks=pretty(1:len, n=10)/len, labels=pretty(1:len, n=10) ) 

			# crest segments
			if ( is.data.frame(crestSub) ) {

				crestSub$start  <- crestSub$start/10/xtotal
				crestSub$end    <- crestSub$end/10/xtotal
			
				crestSub$ymaxcov <- replicate( nrow(crestSub), ymaxcov )
				crestList = split(crestSub, crestSub$type) 

				if (length(crestList$DUP) > 0){
					p <- p + geom_arch( data=crestList$DUP, aes(x=start, xend=end, height=1.2, y=ymaxcov), col='blue', lwd = 0.3)
				}
	          
				if (length(crestList$DEL) > 0){
					p <- p + geom_arch( data=crestList$DEL, aes(x=start, xend=end, height=1.2, y=ymaxcov), col='red', lwd = 0.3)
				}
						 
				if (length(crestList$INV) > 0){
					p <- p + geom_arch( data=crestList$INV, aes(x=start, xend=end, height=1.2, y=ymaxcov), col='purple', lwd = 0.3)
				}

				if (length(crestList$ITX) > 0){ 
					p <- p + geom_arch( data=crestList$ITX, aes( x=start, xend=end, height=1.2, y=ymaxcov ), col="#7CFC00", lwd=0.3 ) 
				}
			  
				if (length(crestList$CTX) > 0){
					p <- p + geom_linerange( data = crestList$CTX, aes(x=start,  y=ymaxcov, ymin=ymaxcov, ymax=ymaxcov+1 ), col="#006400", lty=1, lwd=0.3) 
					for ( i in seq_len( nrow(crestList$CTX) ) ){
						if (crestList$CTX$chromosome[i]==crestList$CTX$chr2[i]){
							p <- p + geom_text(data=crestList$CTX[i,], aes(x=start+0.001, y=ymaxcov+1, label = chr2 ), cex = 2, col = "#006400")
						} else {
							p <- p + geom_text(data=crestList$CTX[i,], aes(x=start+0.001, y=ymaxcov+1, label = chr2), cex = 2, col = "#006400")
						}
					}
				}
			}
		}
		return(p)
}

# plot DH mean values
plotDHmeans <- function( seg, chrLen,  plots='single', p=NULL ) {

		xtotal = chrLen/10
		ymaxDH = 1
		seg$start    	<- (seg$start/10)/xtotal
		seg$end		<- (seg$end/10)/xtotal

		# segment dh means
		p <- ggplot() 
		p <- p + geom_segment(data=seg, aes( x=start, xend=end, y=dhMean, yend=dhMean), size=1, colour="#0000CD", na.rm=TRUE) + theme_bw()

		# labs, scales etc.
		p <- p + xlab('') + ylab('DH segment means') +   xlim( c(0,1) ) + ylim( c(0, ymaxDH ) )
		p <- p + geom_hline(yintercept=c(0, 0.2, 0.4, 0.6, 0.8, 1), col='#C0C0C0', lty='dotted', lwd=0.5)
		p <- p + scale_y_continuous( breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), label=c(0,"", "", "", "", 1) )
		p <- p + theme( axis.title = element_text(size=12), legend.position="none", panel.grid=element_blank() )

		p <- p + scale_x_continuous( breaks=NULL, label=NULL )					#blank

		# segment boundaries,
		if(plots=='single'){
			p <- p + geom_vline( xintercept=seg$start, col='#000000AA', lty=5, lwd=0.2 )
			p <- p + geom_vline( xintercept=seg$end  , col = "#000000DD", lty=5, lwd =0.2 )
		}

		return(p)
}
	



# plot raw BAF values
plotRawBAF <- function(ratio, seg=NULL, chrLen, plots='single', p=NULL){

		xtotal = chrLen/10
		ymaxDH = 1
		len    = chrLen/1000000
		haploColors = c("#000000","red", "blue")
		colScale <- scale_colour_manual( values = c( '0'=haploColors[1], '1'=haploColors[2], '2'=haploColors[3] ) )

		seg$start    	<- (seg$start/10)/xtotal
		seg$end		<- (seg$end/10)/xtotal

		ratio$start	<- (ratio$start/10)/xtotal
		ratio$end	<- (ratio$end/10)/xtotal
		ratio$SNP	<- (ratio$SNP/10)/xtotal

		# BAF values
		p <-  ggplot(environment=environment())
		if (nrow(ratio) >0){
			p <- p + geom_point( data=ratio, aes( SNP, betaT, col=as.character(haplotype)), pch=16, cex=0.3) 
		}
		# labs, scales etc.
		p <- p + geom_hline(yintercept=c(0, 0.2, 0.4, 0.6, 0.8, 1), col="#C0C0C0", lty="dotted", lwd=0.5)
		p <- p + xlab("") + ylab("raw BAF") + ylim(c(0,ymaxDH)) + xlim( c(0,1) ) + theme_bw()# line=6
		p <- p + scale_y_continuous( breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), label=c(0,"","","","",1) )
		p <- p + theme( axis.title = element_text(size=12), legend.position="none", panel.grid=element_blank(), panel.background=element_blank(), panel.margin=unit(c(0,0,0,0), 'mm') )
		p <- p + colScale
		# segment boundaries,
		if ( plots=='single' ){
			p <- p + scale_x_continuous( breaks=pretty(1:len, n=10)/len, labels=pretty(1:len, n=10) )
			p <- p + geom_vline( xintercept=seg$start, col="#000000AA", lty=5, lwd=0.2 )
			p <- p + geom_vline( xintercept=seg$end, col = "#000000DD", lty=5, lwd =0.2 )
		}

		return(p)
}

# adjust values with obtained ploidies and puirities  

completeSNP = function( chr, dat, Ploidy, Purity, fullPloidy ) {
	
	if (nrow(dat) < 1 ){
		GNL <- data.frame( matrix( vector(), 0,1, dimnames=list(c(), 'GNL' ) ), stringsAsFactors=F )
		dat_new <- cbind( dat, GNL )
		return(dat_new)
	}

	D = Purity * Ploidy + 2 * (1 - Purity) 		

	dat$copyT = (dat$copyT * D - 2 * (1-Purity) ) / (Purity)
	dat$meanTCN = (dat$meanTCN * D - 2 * (1-Purity) ) / (Purity)
	dat$roundTCN = round(dat$meanTCN)

	dat$GNL = NA

	if (sex=='male' & (chr==23 | chr == 24)) {
		full_ploi = fullPloidy/2

		dat$copyT = dat$copyT / 2
		dat$meanTCN = dat$meanTCN / 2
		dat$roundTCN = round(dat$meanTCN)
	}else if ( sex=='klinefelter' &  chr == 24 ) {
		full_ploi = fullPloidy/2
		dat$copyT = dat$copyT / 2
		dat$meanTCN = dat$meanTCN / 2
		dat$roundTCN = round(dat$meanTCN)
	}else{
		full_ploi = fullPloidy
	}

	sel = which(dat$roundTCN == full_ploi)
	dat$GNL[sel] = 'n'
	rm(sel)

	sel = which(dat$roundTCN > full_ploi)
	dat$GNL[sel] = 'g'
	rm(sel)

	sel = which(dat$roundTCN < full_ploi)
	dat$GNL[sel] = 'l'
	rm(sel)
	
	return(dat)
}

completeSeg = function( comb, Ploidy, Purity, id, solutionPossible=NA ) {

	comb$tcnMeanRaw <- comb$tcnMean

	#calculate correct TCN with estimated ploidy and purity
	D = Purity * Ploidy + 2 * ( 1 - Purity )
	comb$tcnMean = (comb$tcnMean - ( 2 *( 1-Purity ))/D ) / (Purity/D)

	fullPloidyLength <- sapply(unique(round(comb$tcnMean)), function(i) sum( as.numeric( comb$length[round(comb$tcnMean)==i] ) ) )
	fullPloidyTab <- data.frame( ploidy=unique(round(comb$tcnMean)), length=fullPloidyLength )
	sel <- which(fullPloidyTab$length==max(fullPloidyTab$length))
	fullPloidy <- fullPloidyTab$ploidy[sel[1]]
	if ( length(sel) > 1){
		cat("WARNING two plausible full ploidies found, the selected solution might not be the correct one\n")
	}
	rm(sel)
#	fullPloidy <- which(tabulate(round(comb$tcnMean)) == max(tabulate(round(comb$tcnMean))) )
#	if ( length(fullPloidy) > 1 ){
#		tmp <- fullPloidy - Ploidy
#		sel <- which( tmp == min(tmp) )
#		fullPloidy <- fullPloidy[sel]
#	}

	comb$AF = NA
	comb$BAF = NA
	comb$dhMean = NA
	comb$c1Mean = NA
	comb$c2Mean = NA
	comb$A = NA
	comb$B = NA
	comb$genotype = NA
	comb$TCN = NA
	comb$GNL = NA

	# for balanced alleles split cn equally
	sel = which(comb$peaks == 1 & 
	            comb$meanCovT != "NA" & 
	            comb$tcnMean != "NA" & 
	            comb$meanCovB != "NA" & 
	            comb$meanCovB != "NaN")
	            
	comb$dhMean[sel] = 0
	comb$c1Mean[sel] = comb$tcnMean[sel] / 2
	comb$c2Mean[sel] = comb$tcnMean[sel] / 2
	rm(sel)

	# unbalanced correct with allelic factor
	sel = which(comb$peaks == 2 & 
	            comb$meanCovT != "NA" & 
	            comb$tcnMean != "NA" & 
	            comb$meanCovB != "NA" & 
	            comb$meanCovB != "NaN")

	comb$AF[sel] = (comb$meanCovT[sel]/10000) / ((Purity*comb$tcnMean[sel])+((1-Purity)*2))			#Initial segnments for covT where 10000 bases long
	comb$BAF[sel] = ((comb$meanCovB[sel]/comb$AF[sel])-(1-Purity))/(Purity*comb$tcnMean[sel])
	comb$dhMean[sel] = 2*(abs(comb$BAF[sel] - 0.5))
	sel_dh = which(comb$dhMean > 1)
	comb$dhMean[sel_dh] = 1
	comb$c1Mean[sel] = 0.5*(1-as.numeric(comb$dhMean[sel]))*as.numeric(comb$tcnMean[sel])
	sel_c1 = which(comb$c1Mean < 0)
	comb$c1Mean[sel_c1] = 0
	comb$c2Mean[sel] = as.numeric(comb$tcnMean[sel]) - as.numeric(comb$c1Mean[sel])
	rm(sel)

	sel = which(is.na(comb$dhMean))
	comb$c1Mean[sel] = NA
	comb$c2Mean[sel] = NA
	rm(sel)

	# check whether difference to closest plausible copy number is >0.3
	# if so => assume subpopulation
	sel = which(!is.na(comb$dhMean))
	for ( s in seq_along(sel)){
		if ( abs( round( as.numeric( comb$c2Mean[sel[s]] ),0 ) - as.numeric( comb$c2Mean[sel[s]] ) ) <= 0.3 |
		     abs( round( as.numeric( comb$c2Mean[sel[s]] ),0 ) - as.numeric( comb$c2Mean[sel[s]] ) ) >= 0.7 ){
			comb$A[sel[s]] = round(as.numeric(comb$c2Mean[sel[s]]), 0)
		}else{
			comb$A[sel[s]] = "sub"
		}
		if ( abs( round( as.numeric( comb$c1Mean[sel[s]] ),0 ) - as.numeric( comb$c1Mean[sel[s]] ) ) <= 0.3 |
		     abs( round( as.numeric( comb$c1Mean[sel[s]] ),0 ) - as.numeric( comb$c1Mean[sel[s]] ) ) >= 0.7 ){
			comb$B[sel[s]] = round(as.numeric(comb$c1Mean[sel[s]]), 0)
		}else{
			comb$B[sel[s]] = "sub"
		}
		comb$genotype[sel[s]] = paste(comb$A[sel[s]], comb$B[sel[s]], sep=":")
	}
	rm(sel)           

	sel = which( ! is.na(comb$tcnMean) &
		     ! is.na(comb$A) &
		     ! is.na(comb$B) )

	for (s in seq_along(sel)){
		if ( comb$A[sel[s]] != "sub" & comb$B[sel[s]] != 'sub'){
			comb$TCN[sel[s]] <- as.numeric(comb$A[sel[s]]) + as.numeric(comb$B[sel[s]]) 
		}else {
			comb$TCN[sel[s]] <- "sub"
		}
	}
	rm(sel)

	# classify according to dh and copy numbers
	# all segments with defined peak were heterozygous in the normal tissue
	sel <- which( ! is.na(comb$tcnMean) & 
		      ! is.na(comb$dhMean)  &
		      comb$TCN != "sub" &
		      comb$TCN == 2 &
#		      comb$TCN == fullPloidy &
		      comb$dhMean == 0 )
	comb$GNL[sel] <- 'neutral'
	rm(sel)

	sel <- which( ! is.na(comb$tcnMean) & 
	              ! is.na(comb$dhMean)  &
	      	      comb$TCN != "sub" &
#	      	      comb$TCN == fullPloidy &
		      comb$TCN == 2 &
		      (round(comb$c1Mean) == 0 |
		      round(comb$c2Mean) == 0 ) )

	comb$GNL[sel] <- 'cnLOH'
	rm(sel)

	sel <- which( ! is.na(comb$tcnMean) & 
	              ! is.na(comb$dhMean)  &
	      	      comb$TCN != "sub" &
#	      	      comb$TCN == fullPloidy &
	      	      comb$TCN == 2 &
		    ( round(comb$c1Mean) != 0 &
		      round(comb$c2Mean) != 0   ) & 
		      comb$dhMean > 0 )

	comb$GNL[sel] <- 'mixDH'
	rm(sel)

	sel <- which( ! is.na(comb$tcnMean) & 
	              ! is.na(comb$dhMean)  &
	      	      comb$TCN != "sub" &
#	      	      comb$TCN > fullPloidy &
	      	      comb$TCN > 2 &
		      (round(comb$c1Mean) == 0 |
		      round(comb$c2Mean) == 0 ) )

	comb$GNL[sel] <- 'LOHgain'
	rm(sel)

	sel <- which( ! is.na(comb$tcnMean) & 
	              ! is.na(comb$dhMean)  &
	      	      comb$TCN != "sub" &
#	      	      comb$TCN > fullPloidy &
	      	      comb$TCN > 2 &
		    ( round(comb$c1Mean) != 0 &
		      round(comb$c2Mean) != 0 ) )

	comb$GNL[sel] <- 'gain'
	rm(sel)

	sel <- which( ! is.na(comb$tcnMean) & 
	              ! is.na(comb$dhMean)  &
	      	      comb$TCN != "sub" &
#	      	      comb$TCN < fullPloidy &
	      	      comb$TCN < 2 &
		    ( round(comb$c1Mean) != 0 &
		      round(comb$c2Mean) != 0 ) )

	comb$GNL[sel] <- 'loss'
	rm(sel)

	sel <- which( ! is.na(comb$tcnMean) & 
	              ! is.na(comb$dhMean)  &
	      	      comb$TCN != "sub" &
#	      	      comb$TCN < fullPloidy &
	      	      comb$TCN < 2 &
		      (round(comb$c1Mean) == 0 |
		      round(comb$c2Mean) == 0 ) )

	comb$GNL[sel] <- 'LOH'
	rm(sel)

	sel <- which( comb$TCN == 'sub' )
	comb$GNL[sel] <- 'sub'


	# sex specific
	if (sex == "male") {
		sel = which(comb$chromosome == 23 | comb$chromosome == 24)
		comb$tcnMean[sel] = comb$tcnMean[sel] / 2
		comb$A[sel] = round(comb$tcnMean[sel], 0)
		comb$B[sel] = NA
		comb$AF[sel] = NA
		comb$BAF[sel] = NA
		comb$dhMean[sel] = NA
		comb$c1Mean[sel] = NA
		comb$c2Mean[sel] = NA  
		comb$genotype[sel] = comb$A[sel]
		comb$TCN[sel] = comb$A[sel]
		rm(sel)

#		sel = which((comb$chromosome == 23 | comb$chromosome == 24) & (comb$TCN == (fullPloidy/2)))
		sel = which((comb$chromosome == 23 | comb$chromosome == 24) & (comb$TCN == 1))
		comb$GNL[sel] = "neutral"
		rm(sel)

		sel = which((comb$chromosome == 23 | comb$chromosome == 24) & (comb$TCN > 1 ))
		comb$GNL[sel] = "gain"
		rm(sel)

		sel = which((comb$chromosome == 23 | comb$chromosome == 24) & (comb$TCN < 1))
		comb$GNL[sel] = "loss"
		rm(sel)
	}else if (sex == "klinefelter") {
		sel = which(comb$chromosome == 24)
		comb$tcnMean[sel] = comb$tcnMean[sel] / 2
		comb$A[sel] = round(comb$tcnMean[sel], 0)
		comb$B[sel] = NA
		comb$AF[sel] = NA
		comb$BAF[sel] = NA
		comb$dhMean[sel] = NA
		comb$c1Mean[sel] = NA
		comb$c2Mean[sel] = NA  
		comb$genotype[sel] = comb$A[sel]
		comb$TCN[sel] = comb$A[sel]
		rm(sel)

		sel = which((comb$chromosome == 24) & (comb$TCN == 1 ))
		comb$GNL[sel] = "neutral"
		rm(sel)

		sel = which((comb$chromosome == 24) & (comb$TCN > 1 ))
		comb$GNL[sel] = "gain"
		rm(sel)

		sel = which((comb$chromosome == 24) & (comb$TCN < 1 ))
		comb$GNL[sel] = "loss"
		rm(sel)

	} else if (sex == "female") {
		sel = which(comb$chromosome == 24)
		if (length(sel) >= 1) {
			comb = comb[-sel, ]
		}
		rm(sel)
	}

	sel = which(comb$map == "homozygousDel")
	if (length(sel)>0){
		comb$TCN[sel] = 0
		comb$tcnMean[sel] = 0
		comb$A[sel] = 0
		comb$B[sel] = 0
		comb$dhMean[sel] = 0
		comb$c1Mean[sel] = 0
		comb$c2Mean[sel] = 0
		comb$GNL[sel] = "homozygousDel"
		comb$map[sel] = 'mappable'
	}

	if (length(sel)>0){		  
		sel = which((comb$map == "unmappable") & (comb$chromosome != 23) & (comb$chromosome != 24))
		comb$TCN[sel] = NA
		comb$tcnMean[sel] = NA
		comb$A[sel] = NA
		comb$B[sel] = NA
		comb$dhMean[sel] = NA
		comb$c1Mean[sel] = NA
		comb$c2Mean[sel] = NA
		comb$GNL[sel] = NA
	}

	#format data so that no e-x is used and 0.5 is rounded to the next bigger value for start and next lower for end of a segment
	comb$start	<- as.integer(ceiling(comb$start))
	comb$end	<- as.integer(floor(comb$end))
	comb		<- comb[order(comb$chromosome, comb$start),]
	comb_out   	<- format(comb, scientific = FALSE, trim = TRUE)
	colnames(comb_out)[1] <- "#chromosome"


	write.table(comb_out, qq("@{outDir}/@{id}_comb_pro_extra@{round(Ploidy, digits = 3)}_@{Purity}.txt"), sep = "\t", row.names = FALSE, quote = FALSE) 

	important_cols <- c('#chromosome', 'start', 'end', 'length', 'tcnMeanRaw', 'tcnMean', 'crest', 'c1Mean', 'c2Mean', 'dhMean', 'dhMax', 'genotype', 'GNL', 'tcnNbrOfHets','minStart', 'maxStart', 'minStop', 'maxStop')
  	important_sub  <- comb_out[,important_cols]
  
	colnames(important_sub) <- sub("tcnMeanRaw", "covRatio", colnames(important_sub))
	colnames(important_sub) <- sub("tcnMean", "TCN", colnames(important_sub))
	colnames(important_sub) <- sub("crest", "SV.Type", colnames(important_sub))
	colnames(important_sub) <- sub("dhMean", "dhEst", colnames(important_sub))
	colnames(important_sub) <- sub("dhMax", "dhSNPs", colnames(important_sub))
	colnames(important_sub) <- sub("tcnNbrOfHets", "NbrOfHetsSNPs", colnames(important_sub))

	important_sub 	    <- format(important_sub, scientific = FALSE, trim = TRUE)

	qual = sum( important_sub$GNL != 'sub', na.rm=TRUE) / sum(!is.na(important_sub$GNL))
	importantFile = qq("@{outDir}/@{id}_most_important_info@{round(Ploidy, digits = 3)}_@{Purity}.txt")

	#change parameter names for json conversion output
	normalContamination= 1-purity
	goodnessOfFit=qual
	ploidyFactor=Ploidy
	ploidy=fullPloidy
	caller = "ACEseq"
	gender = sex
	tabFileForJson = qq("@{outDir}/@{id}_cnv_parameter_@{round(Ploidy, digits = 3)}_@{Purity}.txt")
	write.table( data.frame( normalContamination, ploidyFactor, ploidy, goodnessOfFit, gender, solutionPossible ), tabFileForJson, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t' )

	write.table(qq("#purity:@{Purity}\n#ploidy:@{ploidyFactor}\n#roundPloidy:@{fullPloidy}\n#fullPloidy:@{fullPloidy}\n#quality:@{qual}\n#assumed sex:@{sex}"), importantFile, col.names=FALSE, row.names=FALSE, quote=FALSE )
	write.table( important_sub,
		     importantFile,
		     sep = "\t", row.names = FALSE, quote = FALSE, append=TRUE ) 

	return(list(comb,fullPloidy))
}	
