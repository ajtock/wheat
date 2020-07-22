#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Change xblocks height to 46.0 in chrPartitionPlotCov2_feature2 function

# Usage:
# ./chrProfilePlot_1homoeologuegroup_peaks_cMMb_genes_TEsuperfam_step1Mb_altColours.R ASY1_CS ASY1_CS_Rep1_ChIP 10Mb 1Mb darkgreen 200720 'Mariner_DTT' 'chr3A,chr3B,chr3D'

#markChIPA <- "ASY1_CS"
#libNameChIPA <- "ASY1_CS_Rep1_ChIP"
#winName <- "10Mb"
#stepName <- "1Mb"
#colourA <- "darkgreen"
#date <- "200720"
#superfam <- "Mariner_DTT"
#chrName <- unlist(strsplit("chr3A,chr3B,chr3D",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
winName <- args[3]
stepName <- as.numeric(args[4])
colourA <- args[13]
date <- args[15]
superfam <- args[16]
chrName <- unlist(strsplit(args[17],
                           split = ","))

makeTransparent <- function(thisColour, alpha = 210)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}

colourA <- makeTransparent(colourA)

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)
library(plyr)
library(data.table)
library(varhandle)
library(zoo)

plotDir <- "plots/"

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
featureB <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/genes/gene_frequency_per_10Mb_step1Mb.txt"),
                       header = T)
colnames(featureB) <- c("chr", "window", "filt_feature")
featureB <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/TEs/superfamilies/TE_frequency_per_10Mb_step1Mb",
                              "_superfamily_", superfam, ".txt"),
                       header = T)
colnames(featureB) <- c("chr", "window", "filt_feature")


filt_chrProfilescMMb <- mclapply(seq_along(chrs), function(x) {
  cMMb[cMMb$chr == chrs[x],]
}, mc.cores = length(chrs))

#filt_chrProfilescMMb <- mclapply(seq_along(chrProfilescMMb), function(x) {
#  filt_chrProfilecMMb <- stats::filter(x = chrProfilescMMb[[x]]$cMMb,
#                                           filter = f,
#                                           sides = 2)
#  # Given missing cM/Mb data for some of the more distal windows,
#  # need a different way of extending the leftmost and rightmost
#  # non-NA values to the ends of each chromosome, replacing NAs where they are present
#  leftFlank <- which(is.na(filt_chrProfilecMMb))[which(is.na(filt_chrProfilecMMb)) < N*2]
#  rightFlank <- which(is.na(filt_chrProfilecMMb))[which(is.na(filt_chrProfilecMMb)) > N*2]
#  filt_chrProfilecMMb[leftFlank] <- filt_chrProfilecMMb[leftFlank[length(leftFlank)]+1]
#  filt_chrProfilecMMb[rightFlank] <- filt_chrProfilecMMb[rightFlank[1]-1]
##  filt_chrProfilecMMb[1:flank] <- filt_chrProfilecMMb[flank+1]
##  filt_chrProfilecMMb[(length(filt_chrProfilecMMb)-flank+1):length(filt_chrProfilecMMb)] <- filt_chrProfilecMMb[(length(filt_chrProfilecMMb)-flank)]
#  data.frame(chr = as.character(chrProfilescMMb[[x]]$chr),
#             window = as.integer(chrProfilescMMb[[x]]$windowStart),
#             filt_feature = as.numeric(filt_chrProfilecMMb),
#             stringsAsFactors = F)
#}, mc.cores = length(chrProfilescMMb))

filt_chrProfilesFeatureB <- mclapply(seq_along(chrs), function(x) {
  featureB[featureB$chr == chrs[x] &
           featureB[featureB$chr == chrs[x],]$window %in% cMMb[cMMb$chr == chrs[x],]$window,]
}, mc.cores = length(chrs))

#filt_chrProfilesFeatureB <- mclapply(seq_along(chrProfilesFeatureB), function(x) {
#  filt_chrProfileFeatureB <- stats::filter(x = chrProfilesFeatureB[[x]]$cMMb,
#                                           filter = f,
#                                           sides = 2)
#  # Given missing cM/Mb data for some of the more distal windows,
#  # need a different way of extending the leftmost and rightmost
#  # non-NA values to the ends of each chromosome, replacing NAs where they are present
#  leftFlank <- which(is.na(filt_chrProfileFeatureB))[which(is.na(filt_chrProfileFeatureB)) < N*2]
#  rightFlank <- which(is.na(filt_chrProfileFeatureB))[which(is.na(filt_chrProfileFeatureB)) > N*2]
#  filt_chrProfileFeatureB[leftFlank] <- filt_chrProfileFeatureB[leftFlank[length(leftFlank)]+1]
#  filt_chrProfileFeatureB[rightFlank] <- filt_chrProfileFeatureB[rightFlank[1]-1]
##  filt_chrProfileFeatureB[1:flank] <- filt_chrProfileFeatureB[flank+1]
##  filt_chrProfileFeatureB[(length(filt_chrProfileFeatureB)-flank+1):length(filt_chrProfileFeatureB)] <- filt_chrProfileFeatureB[(length(filt_chrProfileFeatureB)-flank)]
#  data.frame(chr = as.character(chrProfilesFeatureB[[x]]$chr),
#             window = as.integer(chrProfilesFeatureB[[x]]$windowStart),
#             filt_feature = as.numeric(filt_chrProfileFeatureB),
#             stringsAsFactors = F)
#}, mc.cores = length(chrProfilesFeatureB))

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
                                xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                dat1A = filt_chrProfilesChIPA[[x]]$filt_feature,
                                col1A = colourA,
                                dat1B = filt_chrProfilescMMb[[x]]$filt_feature,
                                col1B = makeTransparent("darkorange2"),
                                Ylab1 = "",
                                min1A = -max((minCPMA*-1), maxCPMA),
                                max1A = max((minCPMA*-1), maxCPMA),
                                min1B = -max((minCPMB*-1), maxCPMB),
                                max1B = max((minCPMB*-1), maxCPMB),
                                legendLoc = "bottomright",
                                legendLabs = c(sub("_\\w+", "", markChIPA), sub("_\\w+", "", markChIPB), "cM/Mb", "Genes"),
                                xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                col2A = "limegreen",
                                dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                col2B = "magenta3",
                                Ylab2 = "",
                                min2A = 0-maxcMMb,
                                max2A = maxcMMb,
                                min2B = 0-maxFeatureB,
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
                                xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                dat1A = filt_chrProfilesChIPA[[x]]$filt_feature,
                                col1A = colourA,
                                dat1B = filt_chrProfilescMMb[[x]]$filt_feature,
                                col1B = makeTransparent("darkorange2"),
#                                Ylab1 = bquote("Log"[2]*"(ChIP/control)"),
                                Ylab1 = "",
                                min1A = -max((minCPMA_chrs*-1), maxCPMA_chrs),
                                max1A = max((minCPMA_chrs*-1), maxCPMA_chrs),
                                min1B = -max((minCPMB_chrs*-1), maxCPMB_chrs),
                                max1B = max((minCPMB_chrs*-1), maxCPMB_chrs),
                                legendLoc = "bottomright",
                                legendLabs = c(sub("_\\w+", "", markChIPA), "cM/Mb", "Genes", "Mariner TEs"),
                                xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                col2A = "limegreen",
                                dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                col2B = "magenta3",
                                Ylab2 = "",
                                min2A = 0-maxcMMb_chrs,
                                max2A = maxcMMb_chrs,
                                min2B = 0-maxFeatureB_chrs,
                                max2B = maxFeatureB_chrs)
}
dev.off()

