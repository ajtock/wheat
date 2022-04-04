#!/usr/bin/env Rscript

# Compute distance to the centromere of the midpoints of
# physical marker intervals to enable correlations with wild type and fancm mutant
# crossover rates (cM/Mb)

# Usage:
# ./append_distanceToCEN.R

options(stringsAsFactors = F, scipen = 999)
library(data.table)
library(parallel)
library(GenomicRanges)
library(segmentSeq)
library(dplyr)
library(plyr)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

# Table of inter-marker cM values in WT and fancm
#tab <- fread("AxC_markers-mapping-summary_mapped_markers.tsv")
tab <- read.table("AxC_markers-mapping-summary_mapped_markers.tsv",
                  header = T)
tab$pos <- as.integer(gsub(",", "", tab$pos))
tab$chr <- sub("^", "chr", tab$chr)
tab <- tab[ with(tab, base::order(chr, pos)) , ]

interGR <- GRanges()
for(i in unique(tab$chr)) {
  print(i)
  tab_chr <- tab[tab$chr == i,]
  inter_chr_GR <- GRanges(seqnames = c(tab_chr$chr),
                          ranges = IRanges(start = c(1, tab_chr[-nrow(tab_chr),]$pos+1),
                                           end = tab_chr$pos),
                          strand = "*",
                          physical_marker = tab_chr$physical_marker,
                          wt_marker = tab_chr$control_marker,
                          wt_cM = tab_chr$control_cM,
                          fancm_marker = tab_chr$fancm_marker,
                          fancm_cM = tab_chr$fancm_cM) 
  interGR <- c(interGR, inter_chr_GR)
}

#interGR_poswidths <- interGR[width(interGR) > 0]
interGR[width(interGR) == 0]
quantile(width(interGR), probs = seq(0, 1, 0.05))
#          0%           5%          10%          15%          20%          25% 
#        0.00        86.80       263.60       773.20      2502.60      7338.50 
#         30%          35%          40%          45%          50%          55% 
#    22745.90     61824.75    101485.60    145696.75    214491.50    297267.35 
#         60%          65%          70%          75%          80%          85% 
#   413309.40    532912.60    703690.50    970857.50   1340466.00   1869535.45 
#         90%          95%         100% 
#  3301806.00   8655394.25 703285697.00 
#pdf("marker_interval_widths_hist.pdf")
#hist(width(interGR), breaks = 100000, xlim = c(0, 1e6))
#dev.off()


makeDF <- NULL
for(chrName in seqlevels(interGR)) {
  print(chrName)
  chr_interGR <- interGR[seqnames(interGR) == chrName]

  chr_interGR_midpoint <- ( start(chr_interGR)+end(chr_interGR) ) / 2

  chr_CENstart <- centromereStart[which(chrs == chrName)]
  chr_CENend <- centromereEnd[which(chrs == chrName)]

  chr_CENmidpoint <- (chr_CENstart + chr_CENend) / 2

  chr_interGR_distToCEN <- round(abs(chr_interGR_midpoint - chr_CENmidpoint))

  # Make data.frame
  chrDF <- data.frame(chr_interGR,
                      val = chr_interGR_distToCEN)
  # Iteratively combine profiles from each chromosome
  # in a data.frame (profileDNAmeth) that grows with the
  # completion of each loop
  makeDF <- rbind(makeDF, chrDF)
}

makeDF <- makeDF[,c(-5)]
colnames(makeDF) <- c("chr",
                      colnames(makeDF)[2:9],
                      "distToCEN")
makeDF <- makeDF[ with(makeDF, base::order(chr, start, end)) , ]
write.table(makeDF,
            file = paste0("AxC_mapped_marker_intervals_",
                          "distanceToCEN.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
