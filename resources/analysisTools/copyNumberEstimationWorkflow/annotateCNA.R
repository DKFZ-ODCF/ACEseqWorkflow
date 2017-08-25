#!/usr/bin/Rscript

annotateCNA <- function( seg.df, ploidy=fullPloidy, cut.off = 0.7, TCN.colname = "TCN",
                         c1Mean.colname = "c1Mean", c2Mean.colname = "c2Mean", sex=sex ){

       seg.df$CNA.type <- NA
       seg.df$chromosome <- gsub("^23", "X", seg.df$chromosome)
       seg.df$chromosome <- gsub("^24", "Y", seg.df$chromosome)

        selNeutral <- which( ( seg.df[,TCN.colname] >= (ploidy -(1 - cut.off)) &
                                seg.df[,TCN.colname] <= (ploidy + (1 - cut.off)) ) |
                            ( sex=="male" &
                                (seg.df$chromosome =="X" |
                                   seg.df$chromosome== "Y" ) &
                                (seg.df[,TCN.colname] >= (ploidy/2 - (1-cut.off)) &
                                   seg.df[,TCN.colname] <= (ploidy/2 + (1-cut.off)) ) ) |
                            ( sex=="klinefelter" &
                                   seg.df$chromosome== "Y" &
                                (seg.df[,TCN.colname] >= (ploidy/2 - (1-cut.off)) &
                                   seg.df[,TCN.colname] <= (ploidy/2 + (1-cut.off)) ) ) )

        selGain <- which( seg.df[,TCN.colname] >= ploidy + cut.off |
                            ( sex=="male" &
                                (seg.df$chromosome =="X" |
                                   seg.df$chromosome== "Y" ) &
                                seg.df[,TCN.colname] >= (ploidy/2 + cut.off) ) |
                           ( sex=="klinefelter" &
                                   seg.df$chromosome== "Y" &
                                seg.df[,TCN.colname] >= ( ploidy/2 + cut.off) ) )

        selLoss <- which( round(seg.df[,TCN.colname]) > 0 &
                            ( ( seg.df[,TCN.colname] <= ploidy -cut.off &
                                  seg.df$chromosome != "Y"  & 
                                  ( seg.df$chromosome != "X" | sex != "male" ) ) | 
                                ( sex=="male" &
                                    ( seg.df$chromosome == "X" | seg.df$chromosome=="Y") &
                                    seg.df[,TCN.colname] <= (ploidy/2 -cut.off) ) |
                                ( sex=="klinefelter" & seg.df$chromosome=="Y" &
                                    seg.df[,TCN.colname] <= (ploidy/2 -cut.off) ) ) )

        selLOH  <- which( round(seg.df[, c1Mean.colname]) <= 0 & round(seg.df[, c2Mean.colname]) >0 &
                           seg.df$chromosome != "Y" &
                            ( sex =="female" |  sex=="klinefelter" |
                               seg.df$chromosome!="X" ) )

        selHomoDel  <- which( round(seg.df[,TCN.colname]) <= 0 )

        seg.df$CNA.type[selNeutral] <- "TCNneutral"
        seg.df$CNA.type[selGain]    <- paste0( seg.df$CNA.type[selGain], ";DUP")
        seg.df$CNA.type[selLoss]    <- paste0( seg.df$CNA.type[selLoss], ";DEL")
        seg.df$CNA.type[selLOH ]    <- paste0( seg.df$CNA.type[selLOH], ";LOH")
        seg.df$CNA.type[selHomoDel] <- paste0( seg.df$CNA.type[selHomoDel], ";HomoDel")
	return(seg.df)
}
