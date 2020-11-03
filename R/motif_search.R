library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

#' Motif search around PASs
#' 
#' Enrichment of motifs around PASs
#' @param PASs.gr A GRanges object containing exact positions of PASs as single nucleotide region. A metadata column called "PAS.type" is required for each region and should be either sPAS, pPAS or dPAS.
#' @param flanking The number of flanking nucleotides that should be included upstream and downstream.
#' @param motif The motif to be searched.
#' @details For each PAS of a certain PAS type a window is opened in which the desired motif is searched. RNAmaps are generated across all PASs of a certain PAS type.
#' 
#' @return RNAmap plot
#'
#' @import GenomicRanges
#' @import dplyr
#' @import ggplot2
#' @import BSgenome.Mmusculus.UCSC.mm10
#' @import Biostrings
#' @export


motif_search <- function(PASs.gr, flanking=550, motif="GAY"){

  PASs.gr <- PASs.gr + flanking
  
  genome.seq = BSgenome.Mmusculus.UCSC.mm10
  pPAS.motif.start.m <- create.motif.start.matrix(genome.seq, PASs.gr[PASs.gr$new.PAS == "pPAS"], motif = motif, upstream.nt = flanking, downstream.nt = flanking)  
  dPAS.motif.start.m <- create.motif.start.matrix(genome.seq, PASs.gr[PASs.gr$new.PAS == "dPAS"], motif = motif, upstream.nt = flanking, downstream.nt = flanking)
  PAS.motif.start.m <- create.motif.start.matrix(genome.seq, PASs.gr[PASs.gr$new.PAS == "PAS"], motif = motif, upstream.nt = flanking, downstream.nt = flanking)
  
  rna.map.plot.df <- data.frame(x = -550:550,
                                    y = c(colSums(pPAS.motif.start.m) / nrow(pPAS.motif.start.m),
                                          colSums(dPAS.motif.start.m) / nrow(dPAS.motif.start.m),
                                          colSums(PAS.motif.start.m) / nrow(PAS.motif.start.m)),
                                    PAS.type = c(rep("pPAS", ncol(pPAS.motif.start.m)),
                                                 rep("dPAS", ncol(dPAS.motif.start.m)),
                                                 rep("PAS", ncol(PAS.motif.start.m))))
  
  rna.map.plot.df$PAS.type <- factor(rna.map.plot.df$PAS.type, levels = c("pPAS", "dPAS", "PAS"), labels = c("pPAS", "dPAS", "sPAS"))
  
  rna.map.plot <- ggplot(rna.map.plot.df, aes(x=x, y=y, col = PAS.type)) +
    scale_x_continuous(breaks = seq(-500, 500, by = 125)) +
    coord_cartesian(xlim = c(-(flanking-50),(flanking-50))) +
    scale_color_manual(values = c("pPAS" = "#B9529F", "dPAS" = "#009247", "sPAS" = "#58595B")) +
    labs(x = "distance to PAS (nt)", y = "Fraction of PASs") + 
    geom_vline(xintercept = 0, size = 0.25) +
    geom_smooth(size = 1, method = "loess", span=0.1, se = FALSE) +
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



create.motif.start.matrix <- function(genome.seq, PASs.gr, upstream.nt, downstream.nt, motif) {

  PAS.regions.rss <- RNAStringSet(getSeq(genome.seq, PASs.gr))
  
  motif.starts.l <- sapply(vmatchPattern(motif, PAS.regions.rss, fixed=FALSE), start)
  motif.starts.m <- matrix(nrow = length(motif.starts.l), ncol = upstream.nt + downstream.nt +1)
  motif.starts.m[is.na(motif.starts.m)] <- 0
  
  if(length(motif.starts.l) > 0) {
    for (i in 1:length(motif.starts.l)) {
      for (j in motif.starts.l[i]) {
        motif.starts.m[i,j] <- 1
      }
    }
  }
  
  return(motif.starts.m)
}
