#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 01.02.2022

# Examine relationships between differential fancm-wild type recombination rate (cM/Mb) and
# epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq

# Build GLM with differential fancm-wild type cM/Mb as the response variable
# and epigenetic and meiotic protein signals as predictor variables

# Usage:
# ./analyse_differential.R

options(stringsAsFactors = F)
library(fitdistrplus) # descdist, plotdist, fitdist included
library(glm2)
library(MASS) # glm.nb included; MASS is also loaded as a dependency of fitdistrplus
library(pscl) # zeroinfl included
library(vcd) # goodfit included
library(qualityTools) # qqPlot included
library(stats4) # mle (for estimating parameters by maximum likelihood) included
library(dplyr)
#library(segmentSeq)
#library(GenomicRanges)

#library(data.table)
#library(parallel)
#library(GenomicRanges)
#library(dplyr)
#library(plyr)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load cM/Mb data
dat <- read.table("AxC_mapped_marker_intervals_mean_ChIPseq_and_DNAmethyl.tsv", header = T)
colnames(dat)
# [1] "chr"                 "start"               "end"
# [4] "width"               "physical_marker"     "wt_marker"
# [7] "wt_cM"               "fancm_marker"        "fancm_cM"
#[10] "log2_ASY1_CS_input"  "log2_DMC1_input"     "log2_H3K4me1_input"
#[13] "log2_H3K4me3_MNase"  "log2_H3K27ac_input"  "log2_H3K27me3_input"
#[16] "log2_H3K36me3_input" "log2_H3K9me2_MNase"  "log2_H3K27me1_MNase"
#[19] "log2_CENH3_input"    "mCpG"                "mCHG"
#[22] "mCHH"
colnames(dat) <- c(colnames(dat)[1:9],
                   "ASY1", "DMC1", "H3K4me1", "H3K4me3", "H3K27ac",
                   "H3K27me3", "H3K36me3", "H3K9me2", "H3K27me1",
                   "CENH3", "mCG", "mCHG", "mCHH")


dat <- data.frame(dat,
                  


# Inspect distribution of cMMb:
# 1. by plotting the empirical density and the empirical cumulative distribution function (ECDF)
pdf(paste0(plotDir, "cMMb_plotdist.pdf"))
fitdistrplus::plotdist(dat$cMMb, histo = T, demp = T)
dev.off()

# Plot with offset of +1
# "non-positive values not allowed for the 'Gamma' family"
pdf(paste0(plotDir, "cMMb_plus1_plotdist.pdf"))
fitdistrplus::plotdist(dat$cMMb+1, histo = T, demp = T)
dev.off()

# Plot with offset of +1e-06
# "non-positive values not allowed for the 'Gamma' family"
pdf(paste0(plotDir, "cMMb_plus1e-06_plotdist.pdf"))
fitdistrplus::plotdist(dat$cMMb+1e-06, histo = T, demp = T)
dev.off()




