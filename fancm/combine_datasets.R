#!/applications/R/R-4.0.0/bin/Rscript

# Combine signals of epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq,
# within physical marker intervals to enable correlations with wild type and fancm mutant
# crossover rates (cM/Mb)

# Usage:
# ./combine_datasets.R both 1000

#align <- "both"
#genomeBinSize <- 1000

args <- commandArgs(trailingOnly = T)
align <- args[1]
genomeBinSize <- as.integer(args[2])

options(stringsAsFactors = F)
library(data.table)
library(parallel)
library(GenomicRanges)
library(dplyr)
library(plyr)

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
  genomeBinNamePlot <- paste0(genomeBinSize, "-bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e3, "-kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e6, "-Mb")
}

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)


# Load tables of ChIP-seq signal values calculated within marker intervals
ChIP_names <- c(
                "ASY1_CS_Rep1_ChIP_input_SRR6350669",
                "DMC1_Rep1_ChIP_input_SRR6350669",
                "H3K4me1_Rep1_ChIP_SRR8126618_input_SRR6350669",
                "H3K4me3_Rep1_ChIP_MNase_Rep1",
                "H3K27ac_Rep1_ChIP_SRR8126621_input_SRR6350669",
                "H3K27me3_ChIP_SRR6350666_input_SRR6350669",
                "H3K36me3_ChIP_SRR6350670_input_SRR6350669",
                "H3K9me2_Rep1_ChIP_MNase_Rep1",
                "H3K27me1_Rep1_ChIP_MNase_Rep1",
                "CENH3_ChIP_SRR1686799_input_SRR6350669"
               )

ChIP_tab_list <- lapply(seq_along(ChIP_names), function(x) {
  read.table(paste0("AxC_mapped_marker_intervals_", ChIP_names[x], "_",
                    align, "_binSize", genomeBinName, ".tsv"),
             header = T, colClasses = c(rep("NULL", 11), NA))
})

ChIP_DF <- dplyr::bind_cols(ChIP_tab_list)

inter_ChIP_DF <- read.table(paste0("AxC_mapped_marker_intervals_", ChIP_names[1], "_",
                                   align, "_binSize", genomeBinName, ".tsv"),
                           header = T, colClasses = c(rep(NA, 9), rep("NULL", 3)))

# Load tables of context-specific DNA methylation values calculated within marker intervals
DNAmeth_names <- c(
                   "BSseq_Rep8a_SRR6792678_CpG",
                   "BSseq_Rep8a_SRR6792678_CHG",
                   "BSseq_Rep8a_SRR6792678_CHH"
                  )
DNAmeth_tab_list <- lapply(seq_along(DNAmeth_names), function(x) {
  read.table(paste0("AxC_mapped_marker_intervals_", DNAmeth_names[x], ".tsv"),
             header = T, colClasses = c(rep("NULL", 9), NA))
})

DNAmeth_DF <- dplyr::bind_cols(DNAmeth_tab_list)

inter_DNAmeth_DF <- read.table(paste0("AxC_mapped_marker_intervals_", DNAmeth_names[1], ".tsv"),
                               header = T, colClasses = c(rep(NA, 9), rep("NULL", 1)))

stopifnot(identical(inter_ChIP_DF, inter_DNAmeth_DF))
stopifnot(all.equal(inter_ChIP_DF, inter_DNAmeth_DF))

makeDF <- data.frame(inter_ChIP_DF,
                     ChIP_DF,
                     DNAmeth_DF)

write.table(makeDF,
            file = paste0("AxC_mapped_marker_intervals_mean_ChIPseq_and_DNAmethyl.tsv"),
             quote = F, sep = "\t", row.names = F, col.names = T)

