#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Change xblocks ybottom to 213 (depends on superfam ylims) in chrPartitionPlotCovMeth_feature2 function

# Usage:
# ./chrProfilePlot_1homoeologuegroup_DNAmeth_MNase_cMMb_TEsuperfam.R MNase MNase_Rep1 BSseq_Rep8a_SRR6792678 Gypsy_LTR_RLG both 1Mb 1000000 15 darkcyan 120520 'chr3A,chr3B,chr3D'

#markChIPA <- "MNase"
#libNameChIPA <- "MNase_Rep1"
#DNAmethName <- "BSseq_Rep8a_SRR6792678"
#superfam <- "Gypsy_LTR_RLG"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15
#colourA <- "darkcyan"
#date <- "120520"
#chrName <- unlist(strsplit("chr3A,chr3B,chr3D",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
DNAmethName <- args[3]
superfam <- args[4]
align <- args[5]
winName <- args[6]
winSize <- as.numeric(args[7])
N <- as.numeric(args[8])
colourA <- args[9]
date <- args[10]
chrName <- unlist(strsplit(args[11],
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
library(GenomicRanges)

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

DNAmethNamesPlot <- c(
                      "mCG",
                      "mCHG",
                      "mCHH"
                     )
DNAmethColours <- c(
                    "navy",
                    "blue",
                    "deepskyblue1"
                   )
DNAmethColours <- makeTransparent(DNAmethColours)

## ChIPA profile
if(libNameChIPA %in% c("H3K4me3_ChIP_SRR6350668",
                       "H3K27me3_ChIP_SRR6350666",
                       "H3K36me3_ChIP_SRR6350670",
                       "H3K9ac_ChIP_SRR6350667",
                       "CENH3_ChIP_SRR1686799")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else if(libNameChIPA %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                              "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
}
profileChIPA <- read.table(paste0(covDirChIPA, libNameChIPA, "_MappedOn_wheat_v1.0_lowXM_",
                                  align, "_sort_norm_binSize", winName, ".bedgraph"))
if("chrUn" %in% chrs) {
  profileChIPA[dim(profileChIPA)[1],]$V3 <- chrLens[length(chrLens)]
} else {
  profileChIPA[dim(profileChIPA)[1],]$V3 <- 480980714
}
# Rows where the difference between end and start coordinates is > winSize
profileChIPA_bigWins <- profileChIPA[profileChIPA$V3-profileChIPA$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileChIPA <- profileChIPA[profileChIPA$V3-profileChIPA$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileChIPA_bigWinsList <- mclapply(seq_along(1:dim(profileChIPA_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileChIPA_bigWins[x,]$V2,
                      to = profileChIPA_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileChIPA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPA_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileChIPA_bigWins[x,]$V3)),
               V4 = as.numeric(profileChIPA_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileChIPA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPA_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileChIPA_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileChIPA_bigWinsDT <- rbindlist(profileChIPA_bigWinsList)
profileChIPA <- rbind.fill(profileChIPA, profileChIPA_bigWinsDT)
profileChIPA <- profileChIPA[order(profileChIPA$V1, profileChIPA$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPA <- profileChIPA[profileChIPA$V1 == chrs[x],]
  if(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V4))
  }
}, mc.cores = detectCores())
profileChIPA_chrLenValsADT <- rbindlist(chrLenValsAList)
profileChIPA <- rbind.fill(profileChIPA, profileChIPA_chrLenValsADT)
profileChIPA <- profileChIPA[order(profileChIPA$V1, profileChIPA$V2),]

profileChIPA <- data.frame(chr = as.character(profileChIPA$V1),
                           window = as.integer(profileChIPA$V2+1),
                           CPM = as.numeric(profileChIPA$V4),
                           stringsAsFactors = F)

chrProfilesChIPA <- mclapply(seq_along(chrs), function(x) {
  profileChIPA[profileChIPA$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfilesChIPA <- mclapply(seq_along(chrProfilesChIPA), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfilesChIPA[[x]]$CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfilesChIPA[[x]]$chr),
             window = as.integer(chrProfilesChIPA[[x]]$window),
             filt_CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesChIPA))

minCPMA_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))-0.25
maxCPMA_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))+0.25

minCPMA <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))-0.25
maxCPMA <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))+0.25

# DNA methylation proportion profiles
profileDNAmeth <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/DNAmeth/",
                                    DNAmethName, "_per_", winName, ".txt"),
                             header = T)
chrProfilesDNAmeth <- mclapply(seq_along(chrs) , function(x) {
  profileDNAmeth[profileDNAmeth$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfilesDNAmeth <- mclapply(seq_along(chrProfilesDNAmeth), function(x) {
  # mCG
  filt_chrProfile_mCG <- stats::filter(x = chrProfilesDNAmeth[[x]]$mCG,
                                       filter = f,
                                       sides = 2)
  filt_chrProfile_mCG[1:flank] <- filt_chrProfile_mCG[flank+1]
  filt_chrProfile_mCG[(length(filt_chrProfile_mCG)-flank+1):length(filt_chrProfile_mCG)] <- filt_chrProfile_mCG[(length(filt_chrProfile_mCG)-flank)]
  # mCHG
  filt_chrProfile_mCHG <- stats::filter(x = chrProfilesDNAmeth[[x]]$mCHG,
                                        filter = f,
                                        sides = 2)
  filt_chrProfile_mCHG[1:flank] <- filt_chrProfile_mCHG[flank+1]
  filt_chrProfile_mCHG[(length(filt_chrProfile_mCHG)-flank+1):length(filt_chrProfile_mCHG)] <- filt_chrProfile_mCHG[(length(filt_chrProfile_mCHG)-flank)]
  # mCHH
  filt_chrProfile_mCHH <- stats::filter(x = chrProfilesDNAmeth[[x]]$mCHH,
                                        filter = f,
                                        sides = 2)
  filt_chrProfile_mCHH[1:flank] <- filt_chrProfile_mCHH[flank+1]
  filt_chrProfile_mCHH[(length(filt_chrProfile_mCHH)-flank+1):length(filt_chrProfile_mCHH)] <- filt_chrProfile_mCHH[(length(filt_chrProfile_mCHH)-flank)]
  # mC
  filt_chrProfile_mC <- stats::filter(x = chrProfilesDNAmeth[[x]]$mC,
                                      filter = f,
                                      sides = 2)
  filt_chrProfile_mC[1:flank] <- filt_chrProfile_mC[flank+1]
  filt_chrProfile_mC[(length(filt_chrProfile_mC)-flank+1):length(filt_chrProfile_mC)] <- filt_chrProfile_mC[(length(filt_chrProfile_mC)-flank)]
  # Combine in data.frame
  data.frame(chr = as.character(chrProfilesDNAmeth[[x]]$chr),
             window = as.integer(chrProfilesDNAmeth[[x]]$window),
	     filt_mCG = as.numeric(filt_chrProfile_mCG),
             filt_mCHG = as.numeric(filt_chrProfile_mCHG),
             filt_mCHH = as.numeric(filt_chrProfile_mCHH),
             filt_mC = as.numeric(filt_chrProfile_mC),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesDNAmeth))

min_mCG_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesDNAmeth),
  function(x) {
    min(c(filt_chrProfilesDNAmeth[[x]]$filt_mCG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))
max_mCG_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesDNAmeth),
  function(x) {
    max(c(filt_chrProfilesDNAmeth[[x]]$filt_mCG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))

min_mCG <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesDNAmeth[[x]]$filt_mCG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))
max_mCG <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesDNAmeth[[x]]$filt_mCG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))

min_mCHG_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesDNAmeth),
  function(x) {
    min(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))
max_mCHG_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesDNAmeth),
  function(x) {
    max(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))

min_mCHG <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))
max_mCHG <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHG), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))

min_mCHH_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesDNAmeth),
  function(x) {
    min(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHH), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))
max_mCHH_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesDNAmeth),
  function(x) {
    max(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHH), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))

min_mCHH <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHH), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))
max_mCHH <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesDNAmeth[[x]]$filt_mCHH), na.rm = T)
}, mc.cores = length(filt_chrProfilesDNAmeth))))

# Feature frequency chromosome profiles
featureA <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                              "200bp_", winName, ".txt"),
                       header = T)
featureB <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/TEs/superfamilies/TE_frequency_per_",
                              winName, "_superfamily_", superfam, ".txt"),
                       header = T)
colnames(featureB) <- c("chr", "window", "features")

chrProfilesFeatureA <- mclapply(seq_along(chrs), function(x) {
  featureA[featureA$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfilesFeatureA <- mclapply(seq_along(chrProfilesFeatureA), function(x) {
  filt_chrProfileFeatureA <- stats::filter(x = chrProfilesFeatureA[[x]]$cMMb,
                                           filter = f,
                                           sides = 2)
  # Given missing cM/Mb data for some of the more distal windows,
  # need a different way of extending the leftmost and rightmost
  # non-NA values to the ends of each chromosome, replacing NAs where they are present
  leftFlank <- which(is.na(filt_chrProfileFeatureA))[which(is.na(filt_chrProfileFeatureA)) < N*2]
  rightFlank <- which(is.na(filt_chrProfileFeatureA))[which(is.na(filt_chrProfileFeatureA)) > N*2]
  filt_chrProfileFeatureA[leftFlank] <- filt_chrProfileFeatureA[leftFlank[length(leftFlank)]+1]
  filt_chrProfileFeatureA[rightFlank] <- filt_chrProfileFeatureA[rightFlank[1]-1]
#  filt_chrProfileFeatureA[1:flank] <- filt_chrProfileFeatureA[flank+1]
#  filt_chrProfileFeatureA[(length(filt_chrProfileFeatureA)-flank+1):length(filt_chrProfileFeatureA)] <- filt_chrProfileFeatureA[(length(filt_chrProfileFeatureA)-flank)]
  data.frame(chr = as.character(chrProfilesFeatureA[[x]]$chr),
             window = as.integer(chrProfilesFeatureA[[x]]$windowStart),
             filt_feature = as.numeric(filt_chrProfileFeatureA),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesFeatureA))

chrProfilesFeatureB <- mclapply(seq_along(chrs), function(x) {
  featureB[featureB$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfilesFeatureB <- mclapply(seq_along(chrProfilesFeatureB), function(x) {
  filt_chrProfileFeatureB <- stats::filter(x = chrProfilesFeatureB[[x]]$features,
                                           filter = f,
                                           sides = 2)
  filt_chrProfileFeatureB[1:flank] <- filt_chrProfileFeatureB[flank+1]
  filt_chrProfileFeatureB[(length(filt_chrProfileFeatureB)-flank+1):length(filt_chrProfileFeatureB)] <- filt_chrProfileFeatureB[(length(filt_chrProfileFeatureB)-flank)]
  data.frame(chr = as.character(chrProfilesFeatureB[[x]]$chr),
             window = as.integer(chrProfilesFeatureB[[x]]$window),
             filt_feature = as.numeric(filt_chrProfileFeatureB),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesFeatureB))

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
# mCG
pdf(paste0(plotDir, "Wheat_", paste0(chrName, collapse = "_"),
           "_mCG_", DNAmethName,"_", libNameChIPA,
           "_", superfam, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4, width = 10*length(chrName))
par(mfrow = c(1, length(chrName)))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrName)) {
  chrPartitionPlotCovMeth_feature2(chrx = which(chrs %in% chrName),
                                   title = sub("c", "C", chrs[x]),
                                   cenStart = centromereStart[x],
                                   cenEnd = centromereEnd[x],
                                   xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                   dat1A = filt_chrProfilesDNAmeth[[x]]$filt_mCG,
                                   col1A = DNAmethColours[1],
                                   dat1B = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                   col1B = colourA,
                                   Ylab1 = bquote(""),
                                   min1A = min_mCG,
                                   max1A = max_mCG,
                                   min1B = minCPMA,
                                   max1B = maxCPMA,
                                   legendLoc = "topright",
                                   legendLabs = c(
                                                  sub("_\\w+", "", DNAmethNamesPlot[1]),
                                                  sub("_\\w+", "", markChIPA),
                                                  "cM/Mb",
                                                  paste0(sub("_", " ", sub("_\\w\\w\\w$", "", superfam)), "s")),
                                   xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                   dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                   col2A = "darkorange2",
                                   dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                   col2B = "green2",
                                   Ylab2 = "",
                                   min2A = minFeatureA_chrs,
                                   max2A = maxFeatureA_chrs,
                                   min2B = minFeatureB_chrs,
                                   max2B = maxFeatureB_chrs)
}
dev.off()

pdf(paste0(plotDir, "Wheat",
           "_mCG_", DNAmethName,"_", libNameChIPA,
           "_", superfam, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4*7, width = 10*3)
par(mfrow = c(7, 3))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrs)) {
  chrPartitionPlotCovMeth_feature2(chrx = which(chrs %in% chrs),
                                   title = sub("c", "C", chrs[x]),
                                   cenStart = centromereStart[x],
                                   cenEnd = centromereEnd[x],
                                   xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                   dat1A = filt_chrProfilesDNAmeth[[x]]$filt_mCG,
                                   col1A = DNAmethColours[1],
                                   dat1B = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                   col1B = colourA,
                                   Ylab1 = bquote(""),
                                   min1A = min_mCG_chrs,
                                   max1A = max_mCG_chrs,
                                   min1B = minCPMA_chrs,
                                   max1B = maxCPMA_chrs,
                                   legendLoc = "topright",
                                   legendLabs = c(
                                                  sub("_\\w+", "", DNAmethNamesPlot[1]),
                                                  sub("_\\w+", "", markChIPA),
                                                  "cM/Mb",
                                                  paste0(sub("_", " ", sub("_\\w\\w\\w$", "", superfam)), "s")),
                                   xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                   dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                   col2A = "darkorange2",
                                   dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                   col2B = "green2",
                                   Ylab2 = "",
                                   min2A = minFeatureA_chrs,
                                   max2A = maxFeatureA_chrs,
                                   min2B = minFeatureB_chrs,
                                   max2B = maxFeatureB_chrs)
}
dev.off()

# mCHG
pdf(paste0(plotDir, "Wheat_", paste0(chrName, collapse = "_"),
           "_mCHG_", DNAmethName,"_", libNameChIPA,
           "_", superfam, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4, width = 10*length(chrName))
par(mfrow = c(1, length(chrName)))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrName)) {
  chrPartitionPlotCovMeth_feature2(chrx = which(chrs %in% chrName),
                                   title = sub("c", "C", chrs[x]),
                                   cenStart = centromereStart[x],
                                   cenEnd = centromereEnd[x],
                                   xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                   dat1A = filt_chrProfilesDNAmeth[[x]]$filt_mCHG,
                                   col1A = DNAmethColours[2],
                                   dat1B = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                   col1B = colourA,
                                   Ylab1 = bquote(""),
                                   min1A = min_mCHG,
                                   max1A = max_mCHG,
                                   min1B = minCPMA,
                                   max1B = maxCPMA,
                                   legendLoc = "topright",
                                   legendLabs = c(
                                                  sub("_\\w+", "", DNAmethNamesPlot[2]),
                                                  sub("_\\w+", "", markChIPA),
                                                  "cM/Mb",
                                                  paste0(sub("_", " ", sub("_\\w\\w\\w$", "", superfam)), "s")),
                                   xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                   dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                   col2A = "darkorange2",
                                   dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                   col2B = "green2",
                                   Ylab2 = "",
                                   min2A = minFeatureA_chrs,
                                   max2A = maxFeatureA_chrs,
                                   min2B = minFeatureB_chrs,
                                   max2B = maxFeatureB_chrs)
}
dev.off()

pdf(paste0(plotDir, "Wheat",
           "_mCHG_", DNAmethName,"_", libNameChIPA,
           "_", superfam, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4*7, width = 10*3)
par(mfrow = c(7, 3))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrs)) {
  chrPartitionPlotCovMeth_feature2(chrx = which(chrs %in% chrs),
                                   title = sub("c", "C", chrs[x]),
                                   cenStart = centromereStart[x],
                                   cenEnd = centromereEnd[x],
                                   xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                   dat1A = filt_chrProfilesDNAmeth[[x]]$filt_mCHG,
                                   col1A = DNAmethColours[2],
                                   dat1B = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                   col1B = colourA,
                                   Ylab1 = bquote(""),
                                   min1A = min_mCHG_chrs,
                                   max1A = max_mCHG_chrs,
                                   min1B = minCPMA_chrs,
                                   max1B = maxCPMA_chrs,
                                   legendLoc = "topright",
                                   legendLabs = c(
                                                  sub("_\\w+", "", DNAmethNamesPlot[2]),
                                                  sub("_\\w+", "", markChIPA),
                                                  "cM/Mb",
                                                  paste0(sub("_", " ", sub("_\\w\\w\\w$", "", superfam)), "s")),
                                   xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                   dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                   col2A = "darkorange2",
                                   dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                   col2B = "green2",
                                   Ylab2 = "",
                                   min2A = minFeatureA_chrs,
                                   max2A = maxFeatureA_chrs,
                                   min2B = minFeatureB_chrs,
                                   max2B = maxFeatureB_chrs)
}
dev.off()

# mCHH
pdf(paste0(plotDir, "Wheat_", paste0(chrName, collapse = "_"),
           "_mCHH_", DNAmethName,"_", libNameChIPA,
           "_", superfam, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4, width = 10*length(chrName))
par(mfrow = c(1, length(chrName)))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrName)) {
  chrPartitionPlotCovMeth_feature2(chrx = which(chrs %in% chrName),
                                   title = sub("c", "C", chrs[x]),
                                   cenStart = centromereStart[x],
                                   cenEnd = centromereEnd[x],
                                   xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                   dat1A = filt_chrProfilesDNAmeth[[x]]$filt_mCHH,
                                   col1A = DNAmethColours[3],
                                   dat1B = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                   col1B = colourA,
                                   Ylab1 = bquote(""),
                                   min1A = min_mCHH,
                                   max1A = max_mCHH,
                                   min1B = minCPMA,
                                   max1B = maxCPMA,
                                   legendLoc = "topright",
                                   legendLabs = c(
                                                  sub("_\\w+", "", DNAmethNamesPlot[3]),
                                                  sub("_\\w+", "", markChIPA),
                                                  "cM/Mb",
                                                  paste0(sub("_", " ", sub("_\\w\\w\\w$", "", superfam)), "s")),
                                   xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                   dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                   col2A = "darkorange2",
                                   dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                   col2B = "green2",
                                   Ylab2 = "",
                                   min2A = minFeatureA_chrs,
                                   max2A = maxFeatureA_chrs,
                                   min2B = minFeatureB_chrs,
                                   max2B = maxFeatureB_chrs)
}
dev.off()

pdf(paste0(plotDir, "Wheat",
           "_mCHH_", DNAmethName,"_", libNameChIPA,
           "_", superfam, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4*7, width = 10*3)
par(mfrow = c(7, 3))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrs)) {
  chrPartitionPlotCovMeth_feature2(chrx = which(chrs %in% chrs),
                                   title = sub("c", "C", chrs[x]),
                                   cenStart = centromereStart[x],
                                   cenEnd = centromereEnd[x],
                                   xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                   dat1A = filt_chrProfilesDNAmeth[[x]]$filt_mCHH,
                                   col1A = DNAmethColours[3],
                                   dat1B = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                   col1B = colourA,
                                   Ylab1 = bquote(""),
                                   min1A = min_mCHH_chrs,
                                   max1A = max_mCHH_chrs,
                                   min1B = minCPMA_chrs,
                                   max1B = maxCPMA_chrs,
                                   legendLoc = "topright",
                                   legendLabs = c(
                                                  sub("_\\w+", "", DNAmethNamesPlot[3]),
                                                  sub("_\\w+", "", markChIPA),
                                                  "cM/Mb",
                                                  paste0(sub("_", " ", sub("_\\w\\w\\w$", "", superfam)), "s")),
                                   xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                   dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                   col2A = "darkorange2",
                                   dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                   col2B = "green2",
                                   Ylab2 = "",
                                   min2A = minFeatureA_chrs,
                                   max2A = maxFeatureA_chrs,
                                   min2B = minFeatureB_chrs,
                                   max2B = maxFeatureB_chrs)
}
dev.off()
