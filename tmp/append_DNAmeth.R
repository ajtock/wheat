#!/applications/R/R-4.0.0/bin/Rscript

# Compute signals of epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq,
# within physical marker intervals to enable correlations with wild type and fancm mutant
# crossover rates (cM/Mb)

# Usage:
# ./append_DNAmeth.R BSseq_Rep8a_SRR6792678 CpG

#DNAmethName <- "BSseq_Rep8a_SRR6792678"
#context <- "CpG"

args <- commandArgs(trailingOnly = T)
DNAmethName <- args[1]
context <- args[2] 

options(stringsAsFactors = F)
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

if(grepl(pattern = "SRR67926", x = DNAmethName)) {
  DNAmethDir <- "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/BSseq/snakemake_BSseq/coverage/"
} else {
  stop(paste0("DNAmethName is not compatible with the specified coverage path"))
}

# Load processed DNA methylation data
mC <- read.table(gzfile(paste0(DNAmethDir,
                               DNAmethName, "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_",
                               context, ".gz.bismark.cov.gz")),
                 header = F, stringsAsFactors = F,
                 colClasses = c(rep(NA, 2), "NULL", NA, rep("NULL", 2)))
colnames(mC) <- c("chr", "pos", "propMeth")
mC$propMeth <- as.numeric(mC$propMeth/100)

# Function to convert coverage table into GRanges
makeGR <- function(bedgraph) {
  GRanges(seqnames = bedgraph[,1],
          ranges = IRanges(start = bedgraph[,2],
                           end = bedgraph[,2]),
          strand = "*",
          val = bedgraph[,3])
}

mC_GR <- makeGR(bedgraph = mC)

## Function to find coverage windows that overlap marker intervals
#fOverlaps <- function(interGR, datGR) {
#  fOverlaps_obj <- findOverlaps(query = interGR,
#                                subject = datGR,
#                                type = "any",
#                                select = "all",
#                                ignore.strand = T)
#  fOverlaps_obj
#}
#
#fOverlaps_mC <- fOverlaps(interGR = interGR, datGR = mC_GR)

makeDF <- NULL
for(chrName in seqlevels(interGR)) {
  print(chrName)
  chr_interGR <- interGR[seqnames(interGR) == chrName]

  # mC
  chr_mC_GR <- mC_GR[seqnames(mC_GR) == chrName]

#  # Identify overlapping windows and cytosine coordinates in CG context
#  fOverlaps_mC <- fOverlaps(interGR = chr_interGR, datGR = chr_mC_GR)
#
#  # Convert fOverlaps_mC* into list object equivalent to that
#  # generated by segmentSeq::getOverlaps(), in which each
#  # list element corresponds to a genomic window (marker interval) of
#  # sequentially numbered overlapping cytosine coordinates
#  fOverlaps_mC_list <- mclapply(seq_along(unique(queryHits(fOverlaps_mC))),
#                          function(x) {
#                            subjectHits(fOverlaps_mC[queryHits(fOverlaps_mC) == x])
#                          }, mc.cores = detectCores(), mc.preschedule = T)
  # NOTE: findOverlaps approach does not work where one or more query feature range(s)
  # is/are not overlappyed by any subject ranges; these query ranges are not accounted
  # for as list elements in subsequently generated fOverlaps_mC_list list object
  fOverlaps_mC_list <- segmentSeq::getOverlaps(coordinates = chr_interGR,
                                               segments = chr_mC_GR,
                                               overlapType = "overlapping",
                                               whichOverlaps = T,
                                               ignoreStrand = T)
  # Calculate mean methylation proportion in each genomic window (marker interval)
  inter_mC <- sapply(fOverlaps_mC_list, function(x) {
                mean(chr_mC_GR$val[x], na.rm = T)
              })


  # Make data.frame
  chrDF <- data.frame(chr_interGR,
                      val = inter_mC)
  # Iteratively combine profiles from each chromosome
  # in a data.frame (profileDNAmeth) that grows with the
  # completion of each loop
  makeDF <- rbind(makeDF, chrDF)
}

makeDF <- makeDF[,c(-5)]
colnames(makeDF) <- c("chr",
                      colnames(makeDF)[2:9],
                      paste0("m", context))
makeDF <- makeDF[ with(makeDF, base::order(chr, start, end)) , ]
write.table(makeDF,
            file = paste0("AxC_mapped_marker_intervals_",
                          DNAmethName, "_", context, ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
