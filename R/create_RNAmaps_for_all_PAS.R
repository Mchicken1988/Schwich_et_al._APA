library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(rtracklayer)

#' iCLIP RNAmaps for all PAS
#' 
#' Create RNAmaps at sPASs, pPASs and dPASs for two iCLIP libraries.
#' @param PASs.gr A GRanges object containing exact positions of PASs as single nucleotide region. A metadata column called "PAS.type" is required  for each region and should be either sPAS, pPAS or dPAS.
#' @param iCLIP1.plus.bw Path to the BigWig-File of the plus strand for iCLIP library 1.
#' @param iCLIP1.minus.bw Path to the BigWig-File of the minus strand for iCLIP library 1.
#' @param iCLIP2.plus.bw Path to the BigWig-File of the plus strand for iCLIP library 2.
#' @param iCLIP2.minus.bw Path to the BigWig-File of the minus strand for iCLIP library 2.
#' @param upstream Number of upstream nucleotides to include in the RNAmap.
#' @param downstream Number of downstream nucleotides to include in the RNAmap.
#' @details For each PAS type (i.e. sPAS, pPAS and dPAS) RNA maps for two iCLIP libraries are generated in a user-defined window. For comparison of signal differences between the two iCLIP libraries, two proportions Z-tests are performed for each position. Positions with a significant signal difference (adjusted P value <= 0.01) are indicated in black beneath the signals.    
#' 
#' @return RNAmap plot
#'
#' @import GenomicRanges
#' @import dplyr
#' @import ggplot2
#' @import rtracklayer
#' @export


create_RNAmaps_for_all_PAS <- function(PASs.gr, iCLIP1.plus.bw, iCLIP1.minus.bw, iCLIP2.plus.bw, iCLIP2.minus.bw, upstream=450, downstream=150){
  
  ####################
  # Prepare the data #
  ####################
  
  #Open user-defined windows around the PASs 
  ranges(PASs.gr[strand(PASs.gr) == "+"]) <- IRanges(start = start(PASs.gr[strand(PASs.gr) == "+"]) - upstream,
                                                     end = end(PASs.gr[strand(PASs.gr) == "+"]) + downstream)
  
  ranges(PASs.gr[strand(PASs.gr) == "-"]) <- IRanges(start(PASs.gr[strand(PASs.gr) == "-"]) - downstream,
                                                     end = end(PASs.gr[strand(PASs.gr) == "-"]) + upstream)
  
  #Separate plus and minus strand                                                                          
  PASs.plus.gr <- PASs.gr[strand(PASs.gr) == "+"]
  PASs.minus.gr <- PASs.gr[strand(PASs.gr) == "-"]
  
  #Load BigWig-Files as Rles
  iCLIP1.plus.rle <- import(iCLIP1.plus.bw, as="Rle")
  iCLIP1.minus.rle <- import(iCLIP1.minus.bw, as="Rle")
  iCLIP2.plus.rle <- import(iCLIP2.plus.bw, as="Rle")
  iCLIP2.minus.rle <- import(iCLIP2.minus.bw, as="Rle")
  
  
  #################################################
  # Create Crosslink matrices for iCLIP library 1 #
  #################################################
  
  iCLIP1.Xlinks.plus.m <- tryCatch(as.matrix(iCLIP1.plus.rle[PASs.plus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP1.Xlinks.minus.m <- tryCatch(as.matrix(iCLIP1.minus.rle[PASs.minus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP1.Xlinks.minus.m <- iCLIP1.Xlinks.minus.m[,ncol(iCLIP1.Xlinks.minus.m):1]
  iCLIP1.Xlinks.m <- rbind(iCLIP1.Xlinks.plus.m, iCLIP1.Xlinks.minus.m)
  
  #Create Crosslink matrices for iCLIP library 2
  iCLIP2.Xlinks.plus.m <- tryCatch(as.matrix(iCLIP2.plus.rle[PASs.plus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP2.Xlinks.minus.m <- tryCatch(as.matrix(iCLIP2.minus.rle[PASs.minus.gr]), error=function(e) {matrix(, 0, upstream + downstream + 1)})
  iCLIP2.Xlinks.minus.m <- iCLIP2.Xlinks.minus.m[,ncol(iCLIP2.Xlinks.minus.m):1]
  iCLIP2.Xlinks.m <- rbind(iCLIP2.Xlinks.plus.m, iCLIP2.Xlinks.minus.m)
  
  #Set every Crosslink count > 0 to 1
  iCLIP1.Xlinks.m[iCLIP1.Xlinks.m > 0] <- 1
  iCLIP2.Xlinks.m[iCLIP2.Xlinks.m > 0] <- 1 
  
  #Important: Re-order the GRanges object by strand. Now it has the same order like the matrices
  PASs.gr <- c(PASs.plus.gr, PASs.minus.gr)
  
  
  ######################################
  # Perform the two proportions Z-test #
  ######################################
  
  #Create pPAS, dPAS and sPAS matrices for the two-proportions test
  iCLIP1.Xlinks.pPAS.m <- iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),]
  iCLIP1.Xlinks.dPAS.m <- iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),]
  iCLIP1.Xlinks.PAS.m <- iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),]
  iCLIP2.Xlinks.pPAS.m <- iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),]
  iCLIP2.Xlinks.dPAS.m <- iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),]
  iCLIP2.Xlinks.PAS.m <- iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),]
  
  iCLIP1.pPAS.colSums <- colSums(iCLIP1.Xlinks.pPAS.m)
  iCLIP1.dPAS.colSums <- colSums(iCLIP1.Xlinks.dPAS.m)
  iCLIP1.PAS.colSums <- colSums(iCLIP1.Xlinks.PAS.m)
  iCLIP2.pPAS.colSums <- colSums(iCLIP2.Xlinks.pPAS.m)
  iCLIP2.dPAS.colSums <- colSums(iCLIP2.Xlinks.dPAS.m)
  iCLIP2.PAS.colSums <- colSums(iCLIP2.Xlinks.PAS.m)
  
  # Calculate a Scaling factor to normalize the colSums, which accounts for differences in sequencing depth
  sf <- (sum(sum(iCLIP1.plus.rle > 0)) + sum(sum(iCLIP1.minus.rle > 0)))/(sum(sum(iCLIP2.plus.rle > 0)) + sum(sum(iCLIP2.minus.rle > 0)))
  iCLIP1.pPAS.colSums <- iCLIP1.pPAS.colSums/sf
  iCLIP1.dPAS.colSums <- iCLIP1.dPAS.colSums/sf
  iCLIP1.PAS.colSums <- iCLIP1.PAS.colSums/sf
  
  #Perofmr the test and subsequent FDR correcion via BH
  pPAS.adj.pvalues <- sapply(seq_along(iCLIP1.pPAS.colSums), function(i){
    res <- prop.test(x=c(iCLIP1.pPAS.colSums[i], iCLIP2.pPAS.colSums[i]), n=rep(nrow(iCLIP1.Xlinks.pPAS.m),2))
    return(res$p.value)
  }) %>% p.adjust(., method="BH")
  
  dPAS.adj.pvalues <- sapply(seq_along(iCLIP1.dPAS.colSums), function(i){
    res <- prop.test(x=c(iCLIP1.dPAS.colSums[i], iCLIP2.dPAS.colSums[i]), n=rep(nrow(iCLIP1.Xlinks.dPAS.m),2))
    return(res$p.value)
  }) %>% p.adjust(., method="BH")
  
  PAS.adj.pvalues <- sapply(seq_along(iCLIP1.PAS.colSums), function(i){
    res <- prop.test(x=c(iCLIP1.PAS.colSums[i], iCLIP2.PAS.colSums[i]), n=rep(nrow(iCLIP1.Xlinks.PAS.m),2))
    return(res$p.value)
  }) %>% p.adjust(., method="BH")
  
  
  #################
  # Normalization #
  #################
  
  #Normalize by number of nucleotides covered in genome
  iCLIP1.Xlinks.m <- iCLIP1.Xlinks.m / (sum(sum(iCLIP1.plus.rle > 0)) + sum(sum(iCLIP1.minus.rle > 0)))
  iCLIP2.Xlinks.m <- iCLIP2.Xlinks.m / (sum(sum(iCLIP2.plus.rle > 0)) + sum(sum(iCLIP2.minus.rle > 0)))
  
  #Normalize by PAS frequency
  iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),] <- iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),] / sum(PASs.gr$new.PAS == "pPAS")*10^9
  iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),] <- iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),] / sum(PASs.gr$new.PAS == "dPAS")*10^9
  iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),] <- iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),] / sum(PASs.gr$new.PAS == "PAS")*10^9
  
  iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),] <- iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),] / sum(PASs.gr$new.PAS == "pPAS")*10^9
  iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),] <- iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),] / sum(PASs.gr$new.PAS == "dPAS")*10^9
  iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),] <- iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),] / sum(PASs.gr$new.PAS == "PAS")*10^9
  
  
  #################################
  # Dataframes for the final plot #
  #################################
  
  #Generate the plotting dataframe
  rna.map.plot.df <- data.frame(x = -450:150,
                                y = c(colSums(iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),]),
                                      colSums(iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),]),
                                      colSums(iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),]),
                                      colSums(iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),]),
                                      colSums(iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),]),
                                      colSums(iCLIP2.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),])),
                                PAS.type = c(rep("pPAS", ncol(iCLIP1.Xlinks.m)),
                                             rep("dPAS", ncol(iCLIP1.Xlinks.m)),
                                             rep("PAS", ncol(iCLIP1.Xlinks.m)),
                                             rep("pPAS", ncol(iCLIP2.Xlinks.m)),
                                             rep("dPAS", ncol(iCLIP2.Xlinks.m)),
                                             rep("PAS", ncol(iCLIP2.Xlinks.m))),
                                iCLIP.lib = c(rep("iCLIP1", ncol(iCLIP1.Xlinks.m) * 3),
                                              rep("iCLIP2", ncol(iCLIP2.Xlinks.m) * 3)))
  
  rna.map.plot.df$PAS.type <- factor(rna.map.plot.df$PAS.type, levels = c("PAS", "pPAS", "dPAS"), labels =c("sPAS", "pPAS", "dPAS"))
  rna.map.plot.df$iCLIP.lib <- factor(rna.map.plot.df$iCLIP.lib, levels = c("iCLIP1", "iCLIP2"))
  
  #Determine the number of PASs per PAS type 
  n1 <- nrow(iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "pPAS"),])
  n2 <- nrow(iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "dPAS"),])
  n3 <- nrow(iCLIP1.Xlinks.m[which(PASs.gr$new.PAS == "PAS"),])
  
  #Create an annotation dataframe
  anno.text.df <- data.frame(n=c(n1,n2,n3), PAS.type=c("pPAS","dPAS", "PAS"),
                             x=-300, y=10, iCLIP.lib = "iCLIP1")
  anno.text.df$PAS.type <- factor(anno.text.df$PAS.type, levels = c("PAS", "pPAS", "dPAS"), labels =c("sPAS", "pPAS", "dPAS"))
  
  #Create the dataframe with the significant positions
  adj.pvalues.plot <- data.frame(x=-450:150,
                                 PAS.type=c(rep("PAS",length(PAS.adj.pvalues)),
                                            rep("pPAS",length(pPAS.adj.pvalues)),
                                            rep("dPAS",length(dPAS.adj.pvalues))),
                                 adj.pvalue=c(PAS.adj.pvalues, pPAS.adj.pvalues, dPAS.adj.pvalues))
  adj.pvalues.plot <- adj.pvalues.plot[adj.pvalues.plot$adj.pvalue <=0.01,]
  
  adj.pvalues.plot$PAS.type <- factor(adj.pvalues.plot$PAS.type, levels = c("PAS", "pPAS", "dPAS"), labels = c("sPAS", "pPAS", "dPAS"))
  
  ###############
  # RNAmap plot #
  ###############
  
  rna.map.plot <- ggplot(rna.map.plot.df, aes(x=x, y=y, col = iCLIP.lib)) +
    facet_wrap(facets = vars(PAS.type), ncol=3) +
    geom_text(data = anno.text.df, mapping = aes(x = x, y = y, label=paste0("n = ", n)), col = "black") +
    scale_x_continuous(breaks = seq(-400, 100, by = 100)) +
    coord_cartesian(xlim = c(-400,100), ylim = c(0, 10)) +
    scale_color_manual(values = c("iCLIP1" = "#0042FF", "iCLIP2" = "#F79320")) +
    labs(x = "distance to PAS (nt)", y = "Normalized iCLIP signal", title = "Crosslinks (all)") + 
    geom_vline(xintercept = 0, linetype="dashed", size = 0.25) +
    geom_line(alpha=0.6) +
    geom_smooth(size = 1, method = "loess", span=0.1, se = FALSE) +
    geom_point(data = adj.pvalues.plot, aes(x=x, y=0), col = "#4D4D4D", shape = 73, size = 3.5, stroke=0) +
    theme_bw() +
    theme(
      axis.text=element_text(size=7,face="bold"),
      axis.title.x=element_text(size=8,face="bold"),
      axis.title.y=element_text(size=8,face="bold", vjust = 3),
      plot.title = element_text(size=8,face="bold"),
      legend.text = element_text(size=8,face="bold"),
      legend.title = element_text(size=8,face="bold"),
      legend.position = "bottom",
      legend.direction="horizontal")
  
  rna.map.plot
}