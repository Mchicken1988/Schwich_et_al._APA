library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(rtracklayer)

#' iCLIP RNAmaps for regulated PAS
#' 
#' Create RNAmaps at pPASs and dPASs of transcripts with changes in 3'UTR length for two iCLIP libraries.
#' @param PASs.gr A GRanges object containing exact positions of PASs. Required metadata columns are 'PAS.type' as well as 'new.DaPars.SRSF3.regulation', 'new.DaPars.SRSF7.regulation', 'new.DaPars.CPSF6.regulation' and 'new.DaPars.Diff.regulation', which indicate the change in 3'UTR length of the hosting transcript in the different DaPars comparisons. 
#' @param iCLIP1.plus.bw Path to the BigWig-File of the plus strand for iCLIP library 1.
#' @param iCLIP1.minus.bw Path to the BigWig-File of the minus strand for iCLIP library 1.
#' @param iCLIP2.plus.bw Path to the BigWig-File of the plus strand for iCLIP library 2.
#' @param iCLIP2.minus.bw Path to the BigWig-File of the minus strand for iCLIP library 2.
#' @param upstream Number of upstream nucleotides to include in the RNAmap.
#' @param downstream Number of downstream nucleotides to include in the RNAmap.
#' @param DaPars The DaPars comparison to look at. Choices are 'SRSF3' (Srsf3 KD/Ctrl), 'SRSF7' (Srsf7 KD/Ctrl), 'CPSF6' (Cpsf6 KD/Ctrl) and 'Diff' (Diff/Undiff).
#' @param UTR3.length The set of transcripts to look at. Choices are 'shorter' (3'UTRs are getting shorter in the respective DaPars comparison) and 'longer' (3'UTRs are getting longer in the respective DaPars comparison).
#' @details RNAmaps at pPASs and dPASs of transcripts with changes in 3'UTR length are generated for two iCLIP libraries in a user-defined window.
#' In addition, for each PAS type and iCLIP library binding signals at not affected PASs are shown.
#' Significant signal differences (adjusted P value <= 0.01) between regulated and not regulated PASs are indicated in black beneath the signals.     
#' 
#' @return RNAmap plot
#'
#' @import GenomicRanges
#' @import dplyr
#' @import ggplot2
#' @import rtracklayer
#' @export

create_RNAmaps_for_regulated_PAS <- function(PASs.gr, iCLIP1.plus.bw, iCLIP1.minus.bw, iCLIP2.plus.bw, iCLIP2.minus.bw, upstream=450, downstream=150, DaPars=c("SRSF3","SRSF7","CPSF6","Diff"), UTR3.length=c("shorter", "longer")){
  
  ####################
  # Prepare the data #
  ####################
  
  #Extract regulated PASs
  regulated.PASs.gr <- PASs.gr[which(mcols(PASs.gr)[,paste0("new.DaPars.",DaPars,".regulation")] == UTR3.length)]
  
  #Open user-defined windows around the PASs 
  ranges(regulated.PASs.gr[strand(regulated.PASs.gr) == "+"]) <- IRanges(start = start(regulated.PASs.gr[strand(regulated.PASs.gr) == "+"]) - upstream,
                                                     end = end(regulated.PASs.gr[strand(regulated.PASs.gr) == "+"]) + downstream)
  
  ranges(regulated.PASs.gr[strand(regulated.PASs.gr) == "-"]) <- IRanges(start(regulated.PASs.gr[strand(regulated.PASs.gr) == "-"]) - downstream,
                                                     end = end(regulated.PASs.gr[strand(regulated.PASs.gr) == "-"]) + upstream)
  
  #Separate plus and minus strand                                                                          
  regulated.PASs.plus.gr <- regulated.PASs.gr[strand(regulated.PASs.gr) == "+"]
  regulated.PASs.minus.gr <- regulated.PASs.gr[strand(regulated.PASs.gr) == "-"]
  
  #Load BigWig-Files as Rles
  iCLIP1.plus.rle <- import(iCLIP1.plus.bw, as="Rle")
  iCLIP1.minus.rle <- import(iCLIP1.minus.bw, as="Rle")
  iCLIP2.plus.rle <- import(iCLIP2.plus.bw, as="Rle")
  iCLIP2.minus.rle <- import(iCLIP2.minus.bw, as="Rle")
  
  
  #################################################
  # Create Crosslink matrices for iCLIP library 1 #
  #################################################
  
  iCLIP1.Xlinks.plus.m <- tryCatch(as.matrix(iCLIP1.plus.rle[regulated.PASs.plus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP1.Xlinks.minus.m <- tryCatch(as.matrix(iCLIP1.minus.rle[regulated.PASs.minus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP1.Xlinks.minus.m <- iCLIP1.Xlinks.minus.m[,ncol(iCLIP1.Xlinks.minus.m):1]
  iCLIP1.Xlinks.m <- rbind(iCLIP1.Xlinks.plus.m, iCLIP1.Xlinks.minus.m)
  
  #Create Crosslink matrices for iCLIP library 2
  iCLIP2.Xlinks.plus.m <- tryCatch(as.matrix(iCLIP2.plus.rle[regulated.PASs.plus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP2.Xlinks.minus.m <- tryCatch(as.matrix(iCLIP2.minus.rle[regulated.PASs.minus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP2.Xlinks.minus.m <- iCLIP2.Xlinks.minus.m[,ncol(iCLIP2.Xlinks.minus.m):1]
  iCLIP2.Xlinks.m <- rbind(iCLIP2.Xlinks.plus.m, iCLIP2.Xlinks.minus.m)
  
  #Set every Crosslink count > 0 to 1
  iCLIP1.Xlinks.m[iCLIP1.Xlinks.m > 0] <- 1
  iCLIP2.Xlinks.m[iCLIP2.Xlinks.m > 0] <- 1 
  
  #Important: Re-order the GRanges object by strand. Now it has the same order like the matrices
  regulated.PASs.gr <- c(regulated.PASs.plus.gr, regulated.PASs.minus.gr)
  
  #################
  # Normalization #
  #################
  
  #Normalize by number of nucleotides covered in genome
  iCLIP1.Xlinks.m <- iCLIP1.Xlinks.m / (sum(sum(iCLIP1.plus.rle > 0)) + sum(sum(iCLIP1.minus.rle > 0)))
  iCLIP2.Xlinks.m <- iCLIP2.Xlinks.m / (sum(sum(iCLIP2.plus.rle > 0)) + sum(sum(iCLIP2.minus.rle > 0)))
  
  #Normalize by PAS frequency
  iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "pPAS"),] <- iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "pPAS"),] / sum(regulated.PASs.gr$new.PAS == "pPAS")*10^9
  iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "dPAS"),] <- iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "dPAS"),] / sum(regulated.PASs.gr$new.PAS == "dPAS")*10^9
  iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "PAS"),] <- iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "PAS"),] / sum(regulated.PASs.gr$new.PAS == "PAS")*10^9
  
  iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "pPAS"),] <- iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "pPAS"),] / sum(regulated.PASs.gr$new.PAS == "pPAS")*10^9
  iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "dPAS"),] <- iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "dPAS"),] / sum(regulated.PASs.gr$new.PAS == "dPAS")*10^9
  iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "PAS"),] <- iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "PAS"),] / sum(regulated.PASs.gr$new.PAS == "PAS")*10^9
  
  
  #######################
  # Generate background #
  #######################
  
  #Extract non-regulated PASs
  non.regulated.PASs.gr <- PASs.gr[which(!mcols(PASs.gr)[,paste0("new.DaPars.",DaPars,".regulation")] %in% c("shorter", "longer"))]
  
  #Open user-defined windows around the PASs 
  ranges(non.regulated.PASs.gr[strand(non.regulated.PASs.gr) == "+"]) <- IRanges(start = start(non.regulated.PASs.gr[strand(non.regulated.PASs.gr) == "+"]) - upstream,
                                                                                 end = end(non.regulated.PASs.gr[strand(non.regulated.PASs.gr) == "+"]) + downstream)
  
  ranges(non.regulated.PASs.gr[strand(non.regulated.PASs.gr) == "-"]) <- IRanges(start(non.regulated.PASs.gr[strand(non.regulated.PASs.gr) == "-"]) - downstream,
                                                                                 end = end(non.regulated.PASs.gr[strand(non.regulated.PASs.gr) == "-"]) + upstream)
  
  #Separate pPASs and dPASs
  non.regulated.pPASs.gr <- non.regulated.PASs.gr[non.regulated.PASs.gr$new.PAS == "pPAS"]
  non.regulated.dPASs.gr <- non.regulated.PASs.gr[non.regulated.PASs.gr$new.PAS == "dPAS"]
  
  #Randomly sample background, determine iCLIP signals and store results in dataframes
  pPAS.background.df <- create.random.background.dataframe(non.regulated.pPASs.gr, PAS.type = "pPAS", iCLIP1.plus.rle, iCLIP1.minus.rle, iCLIP2.plus.rle, iCLIP2.minus.rle, reps = 5, setsize = sum(regulated.PASs.gr$new.PAS == "pPAS"))
  dPAS.background.df <- create.random.background.dataframe(non.regulated.dPASs.gr, PAS.type = "dPAS", iCLIP1.plus.rle, iCLIP1.minus.rle, iCLIP2.plus.rle, iCLIP2.minus.rle, reps = 5, setsize = sum(regulated.PASs.gr$new.PAS == "dPAS"))
  
  
  
  #################################
  # Dataframe for the final plot #
  #################################
  
  rna.map.plot.df <- data.frame(x = -450:150,
                                y= c(colSums(iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "pPAS"),]),
                                     colSums(iCLIP1.Xlinks.m[which(regulated.PASs.gr$new.PAS == "dPAS"),]),
                                     colSums(iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "pPAS"),]),
                                     colSums(iCLIP2.Xlinks.m[which(regulated.PASs.gr$new.PAS == "dPAS"),])),
                                sd = 0,
                                PAS.type = c(rep("pPAS", ncol(iCLIP1.Xlinks.m)),
                                             rep("dPAS", ncol(iCLIP1.Xlinks.m)),
                                             rep("pPAS", ncol(iCLIP2.Xlinks.m)),
                                             rep("dPAS", ncol(iCLIP2.Xlinks.m))),
                                iCLIP.lib = c(rep("iCLIP1", ncol(iCLIP1.Xlinks.m) * 2),
                                              rep("iCLIP2", ncol(iCLIP1.Xlinks.m) * 2)),
                                iCLIP.type = "real",
                                zscore = 0,
                                pval = 1,
                                adj.pval = 1)
  
  rna.map.plot.df <- rbind(rna.map.plot.df, pPAS.background.df, dPAS.background.df)
  
  
  #########################
  # Z-scores and P values #
  #########################
  
  #Calculate z-scores
  #iCLIP1 pPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$zscore <-
    (rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$y -
       rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "pPAS",]$y) /
    rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "pPAS",]$sd
  
  #iCLIP1 dPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$zscore <-
    (rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$y -
       rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "dPAS",]$y) /
    rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "dPAS",]$sd
  
  #iCLIP2 pPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$zscore <-
    (rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$y -
       rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "pPAS",]$y) /
    rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "pPAS",]$sd
  #iCLIP2 dPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$zscore <-
    (rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$y -
       rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "dPAS",]$y) /
    rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "background" & rna.map.plot.df$PAS.type == "dPAS",]$sd
  
  #Calculate p-values
  #iCLIP1 pPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$pval <-
    2*pnorm(-abs(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$zscore))
  #iCLIP1 dPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$pval <-
    2*pnorm(-abs(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$zscore))
  #iCLIP2 pPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$pval <-
    2*pnorm(-abs(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$zscore))
  #iCLIP2 dPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$pval <-
    2*pnorm(-abs(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$zscore))
  
  #Calculate adjusted p-values
  #iCLIP1 pPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$adj.pval <-
    p.adjust(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$pval, method = "BH")
  #iCLIP1 dPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$adj.pval <-
    p.adjust(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP1" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$pval, method = "BH")
  #iCLIP2 pPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$adj.pval <-
    p.adjust(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "pPAS",]$pval, method = "BH")
  #iCLIP2 dPAS
  rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$adj.pval <-
    p.adjust(rna.map.plot.df[rna.map.plot.df$iCLIP.lib == "iCLIP2" & rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$PAS.type == "dPAS",]$pval, method = "BH")

  #####################################################################
  # Generate the remaining dataframes that are necessary for the plot #    
  #####################################################################

  rna.map.plot.df$PAS.type <- factor(rna.map.plot.df$PAS.type, levels = c("pPAS", "dPAS"))
  rna.map.plot.df$iCLIP.lib <- factor(rna.map.plot.df$iCLIP.lib, levels = c("iCLIP1", "iCLIP2"))
  rna.map.plot.df$iCLIP.type <- factor(rna.map.plot.df$iCLIP.type, levels = c("real", "background"))
  
  ###############
  # RNAmap plot #
  ###############
  
  rna.map.plot <- ggplot(rna.map.plot.df, aes(x=x, y=y)) +
    facet_wrap(facets = vars(iCLIP.lib:PAS.type), ncol = 2, nrow =2) +
    coord_cartesian(xlim = c(-(upstream-50), (downstream-50))) +
    scale_color_manual(values = c("iCLIP1:real" = "#0042FF", "iCLIP2:real" = "#F79320", "iCLIP1:background" = "grey","iCLIP2:background" = "grey")) +
    scale_fill_manual(values = c("iCLIP1:background" = "grey","iCLIP2:background" = "grey")) +
    labs(x = "distance to PAS (nt)", y = "Normalized iCLIP signal") + 
    geom_vline(aes(xintercept = 0), linetype="dashed", size = 0.25) +
    geom_ribbon(data = rna.map.plot.df[rna.map.plot.df$iCLIP.type  == "background",], aes(ymin = y-sd, ymax = y + sd, fill = iCLIP.lib:iCLIP.type), alpha = 0.5) + 
    geom_line(data = rna.map.plot.df[rna.map.plot.df$iCLIP.type  == "background",], aes(col = iCLIP.lib:iCLIP.type)) + 
    geom_line(data = rna.map.plot.df[rna.map.plot.df$iCLIP.type == "real",], alpha=0.6, aes(col = iCLIP.lib:iCLIP.type)) +
    geom_smooth(data = rna.map.plot.df[rna.map.plot.df$iCLIP.type == "real",], aes(col = iCLIP.lib:iCLIP.type), size=1, se = F, span = 0.1) +
    geom_point(data = rna.map.plot.df[rna.map.plot.df$iCLIP.type == "real" & rna.map.plot.df$adj.pval <= 0.01,], aes(x=x, y=-1), col = "#4D4D4D", shape = 73, size = 3.5, stroke=0) +
    theme_bw() +
    theme(
      axis.text=element_text(size=7,face="bold"),
      axis.title.x=element_text(size=8,face="bold"),
      axis.title.y=element_text(size=8,face="bold", vjust = 3),
      plot.title = element_text(size=8,face="bold"),
      legend.text = element_text(size=8,face="bold"),
      legend.title = element_text(size=8,face="bold"),
      legend.position = "bottom",
      legend.direction="vertical")
  
  rna.map.plot
}



create.random.background.dataframe <- function(non.regulated.PASs.gr, PAS.type, iCLIP1.plus.rle, iCLIP1.minus.rle, iCLIP2.plus.rle, iCLIP2.minus.rle, reps, setsize){
  iCLIP1.Xlinks.reps.m <- matrix(,0,unique(width(non.regulated.PASs.gr)))
  iCLIP2.Xlinks.reps.m <- matrix(,0,unique(width(non.regulated.PASs.gr)))
  
  set.seed("123")
  for (i in 1:reps){
    random.indices <- sample(1:length(non.regulated.PASs.gr), setsize)
    random.picks.gr <- non.regulated.PASs.gr[random.indices]
    
    random.picks.plus.gr <- random.picks.gr[strand(random.picks.gr) == "+"]
    random.picks.minus.gr <- random.picks.gr[strand(random.picks.gr) == "-"]
    
    #Create X-link matrix iCLIP1
    iCLIP1.Xlinks.plus.m <- tryCatch(as.matrix(iCLIP1.plus.rle[random.picks.plus.gr]), error=function(e) {matrix(, 0, 450 + 150 + 1)})
    iCLIP1.Xlinks.minus.m <- tryCatch(as.matrix(iCLIP1.minus.rle[random.picks.minus.gr]), error=function(e) {matrix(, 0, 450 + 150 + 1)})
    iCLIP1.Xlinks.minus.m <- iCLIP1.Xlinks.minus.m[,ncol(iCLIP1.Xlinks.minus.m):1]
    iCLIP1.Xlinks.m <- rbind(iCLIP1.Xlinks.plus.m, iCLIP1.Xlinks.minus.m)
    
    #Create X-link matrix iCLIP2
    iCLIP2.Xlinks.plus.m <- tryCatch(as.matrix(iCLIP2.plus.rle[random.picks.plus.gr]), error=function(e) {matrix(, 0, 450 + 150 + 1)})
    iCLIP2.Xlinks.minus.m <- tryCatch(as.matrix(iCLIP2.minus.rle[random.picks.minus.gr]), error=function(e) {matrix(, 0, 450 + 150 + 1)})
    iCLIP2.Xlinks.minus.m <- iCLIP2.Xlinks.minus.m[,ncol(iCLIP2.Xlinks.minus.m):1]
    iCLIP2.Xlinks.m <- rbind(iCLIP2.Xlinks.plus.m, iCLIP2.Xlinks.minus.m)
    
    
    iCLIP1.Xlinks.m[iCLIP1.Xlinks.m > 0] <- 1
    iCLIP2.Xlinks.m[iCLIP2.Xlinks.m > 0] <- 1 
    
    iCLIP1.Xlinks.m <- iCLIP1.Xlinks.m / (sum(sum(iCLIP1.plus.rle > 0)) + sum(sum(iCLIP1.minus.rle > 0)))
    iCLIP2.Xlinks.m <- iCLIP2.Xlinks.m / (sum(sum(iCLIP2.plus.rle > 0)) + sum(sum(iCLIP2.minus.rle > 0)))
    
    #Important: Now the order in the GRanges object is the same like in the Xlinks matrices
    random.picks.gr <- c(random.picks.plus.gr, random.picks.minus.gr)
    
    #Normalize by PAS frequency
    iCLIP1.Xlinks.m <- iCLIP1.Xlinks.m / length(random.picks.gr) * 10^9
    iCLIP2.Xlinks.m <- iCLIP2.Xlinks.m / length(random.picks.gr) * 10^9
    
    #Metaprofile
    iCLIP1.Xlinks.m <- colSums(iCLIP1.Xlinks.m)
    iCLIP2.Xlinks.m <- colSums(iCLIP2.Xlinks.m)
    
    iCLIP1.Xlinks.reps.m <- rbind(iCLIP1.Xlinks.reps.m, iCLIP1.Xlinks.m)
    iCLIP2.Xlinks.reps.m <- rbind(iCLIP2.Xlinks.reps.m, iCLIP2.Xlinks.m)
  }
  
  random.background.df <- data.frame(x = -450:150,
                                     y = c(apply(iCLIP1.Xlinks.reps.m, 2, mean),apply(iCLIP2.Xlinks.reps.m, 2, mean)),
                                     sd = c(apply(iCLIP1.Xlinks.reps.m, 2, sd),apply(iCLIP2.Xlinks.reps.m, 2, sd)),
                                     PAS.type = PAS.type,
                                     iCLIP.lib = c(rep("iCLIP1", length(-450:150)), rep("iCLIP2", length(-450:150))),
                                     iCLIP.type = c(rep("background", length(-450:150)), rep("background", length(-450:150))),
                                     zscore = 0,
                                     pval = 1,
                                     adj.pval = 1)
  
  return(random.background.df)
  
}
