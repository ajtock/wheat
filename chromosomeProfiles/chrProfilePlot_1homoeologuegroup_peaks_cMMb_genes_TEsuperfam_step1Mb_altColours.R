#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Change xblocks height to 10.0 in chrPartitionPlotCov2_feature2 function

# Usage:
# ./chrProfilePlot_1homoeologuegroup_peaks_cMMb_genes_TEsuperfam_step1Mb_altColours.R DMC1 DMC1_Rep1_ChIP 10Mb 1Mb limegreen 200720 'Mariner_DTT' 'chr3A,chr3B,chr3D'

#markChIPA <- "DMC1"
#libNameChIPA <- "DMC1_Rep1_ChIP"
#winName <- "10Mb"
#stepName <- "1Mb"
#colourA <- "limegreen"
#date <- "200720"
#superfam <- "Mariner_DTT"
#chrName <- unlist(strsplit("chr3A,chr3B,chr3D",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
winName <- args[3]
stepName <- args[4]
colourA <- args[5]
date <- args[6]
superfam <- args[7]
chrName <- unlist(strsplit(args[8],
                           split = ","))

makeTransparent <- function(thisColour, alpha = 210)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}

#colourA <- makeTransparent(colourA)

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)
library(plyr)
library(data.table)
library(varhandle)
library(zoo)

plotDir <- "plots/Figure2/v200720/"

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
markers <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_mapping_data.txt",
                      header = TRUE)

# Peak frequency chromosome profiles
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/peaks/",
                           libNameChIPA, "_peak_frequency_per_",
                           winName, "_step", stepName, ".txt"),
                    header = T) 
colnames(peaks) <- c("chr", "window", "filt_feature")
# Feature frequency chromosome profiles
#cMMb <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
#                              "200bp_10Mb_step1Mb.txt"),
#                       header = T)
cMMb <- read.table("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/iwgsc_refseqv1.0_recombination_rate.txt",
                   header = T)
colnames(cMMb) <- c("chr", "window", "intervalEnd", "nbOfSnps", "filt_feature")
#cMMb <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/cMMb_WGIN_CSxParagon_mapping_data_minInterMarkerDist",
#                              "200bp_10Mb_step1Mb_alt.txt"),
#                       header = T)
#colnames(cMMb) <- c("chr", "window", "winEnd", "filt_feature")
featureA <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/genes/gene_frequency_per_10Mb_step1Mb.txt"),
                       header = T)
colnames(featureA) <- c("chr", "window", "filt_feature")
featureB <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/TEs/superfamilies/TE_frequency_per_10Mb_step1Mb",
                              "_superfamily_", superfam, ".txt"),
                       header = T)
colnames(featureB) <- c("chr", "window", "filt_feature")

filt_chrProfilescMMb <- mclapply(seq_along(chrs), function(x) {
  cMMb[cMMb$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfilesPeaks <- mclapply(seq_along(chrs), function(x) {
  peaks[peaks$chr == chrs[x] &
        peaks$window %in% cMMb[cMMb$chr == chrs[x],]$window,]
}, mc.cores = length(chrs))

filt_chrProfilesFeatureA <- mclapply(seq_along(chrs), function(x) {
  featureA[featureA$chr == chrs[x] &
           featureA$window %in% cMMb[cMMb$chr == chrs[x],]$window,]
}, mc.cores = length(chrs))

filt_chrProfilesFeatureB <- mclapply(seq_along(chrs), function(x) {
  featureB[featureB$chr == chrs[x] &
           featureB$window %in% cMMb[cMMb$chr == chrs[x],]$window,]
}, mc.cores = length(chrs))

minPeaks_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesPeaks),
  function(x) {
    min(c(filt_chrProfilesPeaks[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesPeaks))))
maxPeaks_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesPeaks),
  function(x) {
    max(c(filt_chrProfilesPeaks[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesPeaks))))

minPeaks <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesPeaks[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesPeaks))))
maxPeaks <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesPeaks[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesPeaks))))

mincMMb_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilescMMb),
  function(x) {
    min(c(filt_chrProfilescMMb[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilescMMb))))
maxcMMb_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilescMMb),
  function(x) {
    max(c(filt_chrProfilescMMb[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilescMMb))))

mincMMb <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilescMMb[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilescMMb))))
maxcMMb <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilescMMb[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilescMMb))))

minFeatureA_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesFeatureA),
  function(x) {
    min(c(filt_chrProfilesFeatureA[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureA))))
maxFeatureA_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesFeatureA),
  function(x) {
    max(c(filt_chrProfilesFeatureA[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureA))))

minFeatureA <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesFeatureA[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureA))))
maxFeatureA <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesFeatureA[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureA))))

minFeatureB_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesFeatureB),
  function(x) {
    min(c(filt_chrProfilesFeatureB[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureB))))
maxFeatureB_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesFeatureB),
  function(x) {
    max(c(filt_chrProfilesFeatureB[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureB))))

minFeatureB <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesFeatureB[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureB))))
maxFeatureB <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesFeatureB[[x]]$filt_feature), na.rm = T)
}, mc.cores = length(filt_chrProfilesFeatureB))))

# Plot
pdf(paste0(plotDir, "Wheat_", paste0(chrName, collapse = "_"),
           "_", libNameChIPA, "_peaks_",
           "_featureFreq_chrPlot_winSize", winName, "_step", stepName,
           "_CSxRenan_IWGSCanalysis_v", date, ".pdf"),
    height = 4, width = 10*length(chrName))
par(mfrow = c(1, length(chrName)))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrName)) {
  chrPartitionPlotCov2_feature2(chrx = which(chrs %in% chrName),
                                title = sub("c", "C", chrs[x]),
                                cenStart = centromereStart[x],
                                cenEnd = centromereEnd[x],
                                rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
                                rug1Col = "grey40",
                                xplot1 = filt_chrProfilesPeaks[[x]]$window,
                                dat1A = filt_chrProfilesPeaks[[x]]$filt_feature,
                                col1A = colourA,
                                dat1B = filt_chrProfilescMMb[[x]]$filt_feature,
                                col1B = "darkorange2",
                                Ylab1 = "",
                                min1A = 0,
                                max1A = maxPeaks,
                                min1B = 0,
                                max1B = maxcMMb,
                                legendLoc = "topright",
                                legendLabs = c(paste0(sub("_\\w+", "", markChIPA), " peaks"), "cM/Mb", "Genes", "Mariner TEs"),
                                xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                col2A = "magenta3",
                                dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                col2B = "navy",
                                Ylab2 = "",
                                min2A = 0,
                                max2A = maxFeatureA,
                                min2B = 0,
                                max2B = maxFeatureB)
}
dev.off()

pdf(paste0(plotDir, "Wheat",
           "_", libNameChIPA, "_peaks_",
           "_featureFreq_chrPlot_winSize", winName, "_step", stepName,
           "_CSxRenan_IWGSCanalysis_v", date, ".pdf"),
    height = 4*7, width = 8*3)
par(mfrow = c(7, 3))
par(mar = c(5.0, 6.0, 2.1, 6.0))
for(x in which(chrs %in% chrs)) {
  chrPartitionPlotCov2_feature2(chrx = which(chrs %in% chrs),
                                title = sub("c", "C", chrs[x]),
                                cenStart = centromereStart[x],
                                cenEnd = centromereEnd[x],
                                rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
                                rug1Col = "grey40",
                                xplot1 = filt_chrProfilesPeaks[[x]]$window,
                                dat1A = filt_chrProfilesPeaks[[x]]$filt_feature,
                                col1A = colourA,
                                dat1B = filt_chrProfilescMMb[[x]]$filt_feature,
                                col1B = "darkorange2",
                                Ylab1 = "",
                                min1A = 0,
                                max1A = maxPeaks_chrs,
                                min1B = 0,
                                max1B = maxcMMb_chrs,
                                legendLoc = "topright",
                                legendLabs = c(paste0(sub("_\\w+", "", markChIPA), " peaks"), "cM/Mb", "Genes", "Mariner TEs"),
                                xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                col2A = "magenta3",
                                dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                col2B = "navy",
                                Ylab2 = "",
                                min2A = 0,
                                max2A = maxFeatureA_chrs,
                                min2B = 0,
                                max2B = maxFeatureB_chrs)
}
dev.off()
