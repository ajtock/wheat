#!/applications/R/R-4.0.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Change xblocks height to 46.0 in chrPartitionPlotCov2_feature2 function

# Usage:
# ./chrProfilePlot_1homoeologuegroup_log2_ChIP_input_cMMb_genes_mapParamVar_step1Mb_altColours.R DMC1 DMC1_Rep1_ChIP input H3_input_SRR6350669 1Mb 1000000 15 120421 'chr3A,chr3B,chr3D'

#markChIP <- "DMC1"
#libNameChIP <- "DMC1_Rep1_ChIP"
#markControl <- "input"
#libNameControl <- "H3_input_SRR6350669"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15
#date <- "120421"
#chrName <- unlist(strsplit("chr3A,chr3B,chr3D",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
markChIP <- args[1]
libNameChIP <- args[2]
markControl <- args[3]
libNameControl <- args[4]
winName <- args[5]
winSize <- as.numeric(args[6])
N <- as.numeric(args[7])
date <- args[8]
chrName <- unlist(strsplit(args[9],
                           split = ","))

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)
library(plyr)
library(data.table)
library(varhandle)
library(zoo)

makeTransparent <- function(thisColour, alpha = 210)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}

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

## ChIPA profile
covDirChIPA <- paste0("/home/ajt200/analysis/wheat/",
                      markChIP, "/snakemake_ChIPseq/mapped/",
                      "unique/bg/")
profileChIPA <- read.table(paste0(covDirChIPA, libNameChIP, "_MappedOn_wheat_v1.0_lowXM_",
                                  "unique_sort_norm_binSize", winName, ".bedgraph"))
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

## ControlA profile
covDirControlA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                         markControl, "/snakemake_ChIPseq/mapped/unique/bg/")
profileControlA <- read.table(paste0(covDirControlA, libNameControl, "_MappedOn_wheat_v1.0_lowXM_",
                                     "unique_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileControlA_bigWins <- profileControlA[profileControlA$V3-profileControlA$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileControlA <- profileControlA[profileControlA$V3-profileControlA$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileControlA_bigWinsList <- mclapply(seq_along(1:dim(profileControlA_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileControlA_bigWins[x,]$V2,
                      to = profileControlA_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileControlA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlA_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileControlA_bigWins[x,]$V3)),
               V4 = as.numeric(profileControlA_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileControlA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlA_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileControlA_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileControlA_bigWinsDT <- rbindlist(profileControlA_bigWinsList)
profileControlA <- rbind.fill(profileControlA, profileControlA_bigWinsDT)
profileControlA <- profileControlA[order(profileControlA$V1, profileControlA$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileControl <- profileControlA[profileControlA$V1 == chrs[x],]
  if(chrProfileControl[dim(chrProfileControl)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileControl[dim(chrProfileControl)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileControl[dim(chrProfileControl)[1],]$V4))
  }
}, mc.cores = detectCores())
profileControlA_chrLenValsADT <- rbindlist(chrLenValsAList)
profileControlA <- rbind.fill(profileControlA, profileControlA_chrLenValsADT)
profileControlA <- profileControlA[order(profileControlA$V1, profileControlA$V2),]

profileControlA <- data.frame(chr = as.character(profileControlA$V1),
                              window = as.integer(profileControlA$V2+1),
                              CPM = as.numeric(profileControlA$V4),
                              stringsAsFactors = F)

## ChIPA
# Get coverage within each window
profileChIPA <- data.frame(chr = as.character(profileChIPA$chr),
                           window = as.numeric(profileChIPA$window),
                           CPM = as.numeric(log2((profileChIPA$CPM+1)/(profileControlA$CPM+1))),
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


## ChIPB profile
covDirChIPB <- paste0("/home/ajt200/analysis/wheat/",
                      markChIP, "/snakemake_ChIPseq/mapped/",
                      "both/bg/")
profileChIPB <- read.table(paste0(covDirChIPB, libNameChIP, "_MappedOn_wheat_v1.0_lowXM_",
                                  "both_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileChIPB_bigWins <- profileChIPB[profileChIPB$V3-profileChIPB$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileChIPB <- profileChIPB[profileChIPB$V3-profileChIPB$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileChIPB_bigWinsList <- mclapply(seq_along(1:dim(profileChIPB_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileChIPB_bigWins[x,]$V2,
                      to = profileChIPB_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileChIPB_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPB_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileChIPB_bigWins[x,]$V3)),
               V4 = as.numeric(profileChIPB_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileChIPB_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPB_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileChIPB_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileChIPB_bigWinsDT <- rbindlist(profileChIPB_bigWinsList)
profileChIPB <- rbind.fill(profileChIPB, profileChIPB_bigWinsDT)
profileChIPB <- profileChIPB[order(profileChIPB$V1, profileChIPB$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPB <- profileChIPB[profileChIPB$V1 == chrs[x],]
  if(chrProfileChIPB[dim(chrProfileChIPB)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPB[dim(chrProfileChIPB)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPB[dim(chrProfileChIPB)[1],]$V4))
  }
}, mc.cores = detectCores())
profileChIPB_chrLenValsADT <- rbindlist(chrLenValsAList)
profileChIPB <- rbind.fill(profileChIPB, profileChIPB_chrLenValsADT)
profileChIPB <- profileChIPB[order(profileChIPB$V1, profileChIPB$V2),]

profileChIPB <- data.frame(chr = as.character(profileChIPB$V1),
                           window = as.integer(profileChIPB$V2+1),
                           CPM = as.numeric(profileChIPB$V4),
                           stringsAsFactors = F)

## ControlB profile
covDirControlB <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                         markControl, "/snakemake_ChIPseq/mapped/both/bg/")
profileControlB <- read.table(paste0(covDirControlB, libNameControl, "_MappedOn_wheat_v1.0_lowXM_",
                                     "both_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileControlB_bigWins <- profileControlB[profileControlB$V3-profileControlB$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileControlB <- profileControlB[profileControlB$V3-profileControlB$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileControlB_bigWinsList <- mclapply(seq_along(1:dim(profileControlB_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileControlB_bigWins[x,]$V2,
                      to = profileControlB_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileControlB_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlB_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileControlB_bigWins[x,]$V3)),
               V4 = as.numeric(profileControlB_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileControlB_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlB_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileControlB_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileControlB_bigWinsDT <- rbindlist(profileControlB_bigWinsList)
profileControlB <- rbind.fill(profileControlB, profileControlB_bigWinsDT)
profileControlB <- profileControlB[order(profileControlB$V1, profileControlB$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileControl <- profileControlB[profileControlB$V1 == chrs[x],]
  if(chrProfileControl[dim(chrProfileControl)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileControl[dim(chrProfileControl)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileControl[dim(chrProfileControl)[1],]$V4))
  }
}, mc.cores = detectCores())
profileControlB_chrLenValsADT <- rbindlist(chrLenValsAList)
profileControlB <- rbind.fill(profileControlB, profileControlB_chrLenValsADT)
profileControlB <- profileControlB[order(profileControlB$V1, profileControlB$V2),]

profileControlB <- data.frame(chr = as.character(profileControlB$V1),
                              window = as.integer(profileControlB$V2+1),
                              CPM = as.numeric(profileControlB$V4),
                              stringsAsFactors = F)

## ChIPB
# Get coverage within each window
profileChIPB <- data.frame(chr = as.character(profileChIPB$chr),
                           window = as.numeric(profileChIPB$window),
                           CPM = as.numeric(log2((profileChIPB$CPM+1)/(profileControlB$CPM+1))),
                           stringsAsFactors = F)
                        
chrProfilesChIPB <- mclapply(seq_along(chrs), function(x) {
  profileChIPB[profileChIPB$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfilesChIPB <- mclapply(seq_along(chrProfilesChIPB), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfilesChIPB[[x]]$CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfilesChIPB[[x]]$chr),
             window = as.integer(chrProfilesChIPB[[x]]$window),
             filt_CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesChIPB))


## ChIPC profile
covDirChIPC <- paste0("/home/ajt200/analysis/wheat/",
                      markChIP, "/snakemake_ChIPseq_MAPQ0_XM6/mapped/",
                      "both/bg/")
profileChIPC <- read.table(paste0(covDirChIPC, libNameChIP, "_MappedOn_wheat_v1.0_lowXM_",
                                  "both_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileChIPC_bigWins <- profileChIPC[profileChIPC$V3-profileChIPC$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileChIPC <- profileChIPC[profileChIPC$V3-profileChIPC$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileChIPC_bigWinsList <- mclapply(seq_along(1:dim(profileChIPC_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileChIPC_bigWins[x,]$V2,
                      to = profileChIPC_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileChIPC_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPC_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileChIPC_bigWins[x,]$V3)),
               V4 = as.numeric(profileChIPC_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileChIPC_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPC_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileChIPC_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileChIPC_bigWinsDT <- rbindlist(profileChIPC_bigWinsList)
profileChIPC <- rbind.fill(profileChIPC, profileChIPC_bigWinsDT)
profileChIPC <- profileChIPC[order(profileChIPC$V1, profileChIPC$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPC <- profileChIPC[profileChIPC$V1 == chrs[x],]
  if(chrProfileChIPC[dim(chrProfileChIPC)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPC[dim(chrProfileChIPC)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPC[dim(chrProfileChIPC)[1],]$V4))
  }
}, mc.cores = detectCores())
profileChIPC_chrLenValsADT <- rbindlist(chrLenValsAList)
profileChIPC <- rbind.fill(profileChIPC, profileChIPC_chrLenValsADT)
profileChIPC <- profileChIPC[order(profileChIPC$V1, profileChIPC$V2),]

profileChIPC <- data.frame(chr = as.character(profileChIPC$V1),
                           window = as.integer(profileChIPC$V2+1),
                           CPM = as.numeric(profileChIPC$V4),
                           stringsAsFactors = F)

## ControlC profile
covDirControlC <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                         markControl, "/snakemake_ChIPseq_MAPQ0_XM6/mapped/both/bg/")
profileControlC <- read.table(paste0(covDirControlC, libNameControl, "_MappedOn_wheat_v1.0_lowXM_",
                                     "both_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileControlC_bigWins <- profileControlC[profileControlC$V3-profileControlC$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileControlC <- profileControlC[profileControlC$V3-profileControlC$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileControlC_bigWinsList <- mclapply(seq_along(1:dim(profileControlC_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileControlC_bigWins[x,]$V2,
                      to = profileControlC_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileControlC_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlC_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileControlC_bigWins[x,]$V3)),
               V4 = as.numeric(profileControlC_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileControlC_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlC_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileControlC_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileControlC_bigWinsDT <- rbindlist(profileControlC_bigWinsList)
profileControlC <- rbind.fill(profileControlC, profileControlC_bigWinsDT)
profileControlC <- profileControlC[order(profileControlC$V1, profileControlC$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileControl <- profileControlC[profileControlC$V1 == chrs[x],]
  if(chrProfileControl[dim(chrProfileControl)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileControl[dim(chrProfileControl)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileControl[dim(chrProfileControl)[1],]$V4))
  }
}, mc.cores = detectCores())
profileControlC_chrLenValsADT <- rbindlist(chrLenValsAList)
profileControlC <- rbind.fill(profileControlC, profileControlC_chrLenValsADT)
profileControlC <- profileControlC[order(profileControlC$V1, profileControlC$V2),]

profileControlC <- data.frame(chr = as.character(profileControlC$V1),
                              window = as.integer(profileControlC$V2+1),
                              CPM = as.numeric(profileControlC$V4),
                              stringsAsFactors = F)

## ChIPC
# Get coverage within each window
profileChIPC <- data.frame(chr = as.character(profileChIPC$chr),
                           window = as.numeric(profileChIPC$window),
                           CPM = as.numeric(log2((profileChIPC$CPM+1)/(profileControlC$CPM+1))),
                           stringsAsFactors = F)
                        
chrProfilesChIPC <- mclapply(seq_along(chrs), function(x) {
  profileChIPC[profileChIPC$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfilesChIPC <- mclapply(seq_along(chrProfilesChIPC), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfilesChIPC[[x]]$CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfilesChIPC[[x]]$chr),
             window = as.integer(chrProfilesChIPC[[x]]$window),
             filt_CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesChIPC))



# Define plot ylims
minCPM_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))-0.05
#minCPM_chrs <- -0.4
maxCPM_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))+0.05
#maxCPM_chrs <- 0.4

minCPM <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))-0.05
#minCPM <- -0.4
maxCPM <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))+0.05
#maxCPM <- 0.4

# Feature frequency chromosome profiles
#featureA <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
#                              "200bp_10Mb_step1Mb.txt"),
#                       header = T)
featureA <- read.table("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/iwgsc_refseqv1.0_recombination_rate.txt",
                       header = T)
colnames(featureA) <- c("chr", "window", "intervalEnd", "nbOfSnps", "filt_feature")
#featureA <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/cMMb_WGIN_CSxParagon_mapping_data_minInterMarkerDist",
#                              "200bp_10Mb_step1Mb_alt.txt"),
#                       header = T)
#colnames(featureA) <- c("chr", "window", "winEnd", "filt_feature")
featureB <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/genes/gene_frequency_per_10Mb_step1Mb.txt"),
                       header = T)
colnames(featureB) <- c("chr", "window", "filt_feature")

filt_chrProfilesFeatureA <- mclapply(seq_along(chrs), function(x) {
  featureA[featureA$chr == chrs[x],]
}, mc.cores = length(chrs))

#filt_chrProfilesFeatureA <- mclapply(seq_along(chrProfilesFeatureA), function(x) {
#  filt_chrProfileFeatureA <- stats::filter(x = chrProfilesFeatureA[[x]]$cMMb,
#                                           filter = f,
#                                           sides = 2)
#  # Given missing cM/Mb data for some of the more distal windows,
#  # need a different way of extending the leftmost and rightmost
#  # non-NA values to the ends of each chromosome, replacing NAs where they are present
#  leftFlank <- which(is.na(filt_chrProfileFeatureA))[which(is.na(filt_chrProfileFeatureA)) < N*2]
#  rightFlank <- which(is.na(filt_chrProfileFeatureA))[which(is.na(filt_chrProfileFeatureA)) > N*2]
#  filt_chrProfileFeatureA[leftFlank] <- filt_chrProfileFeatureA[leftFlank[length(leftFlank)]+1]
#  filt_chrProfileFeatureA[rightFlank] <- filt_chrProfileFeatureA[rightFlank[1]-1]
##  filt_chrProfileFeatureA[1:flank] <- filt_chrProfileFeatureA[flank+1]
##  filt_chrProfileFeatureA[(length(filt_chrProfileFeatureA)-flank+1):length(filt_chrProfileFeatureA)] <- filt_chrProfileFeatureA[(length(filt_chrProfileFeatureA)-flank)]
#  data.frame(chr = as.character(chrProfilesFeatureA[[x]]$chr),
#             window = as.integer(chrProfilesFeatureA[[x]]$windowStart),
#             filt_feature = as.numeric(filt_chrProfileFeatureA),
#             stringsAsFactors = F)
#}, mc.cores = length(chrProfilesFeatureA))

filt_chrProfilesFeatureB <- mclapply(seq_along(chrs), function(x) {
  featureB[featureB$chr == chrs[x] &
           featureB$window %in% featureA[featureA$chr == chrs[x],]$window,]
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
           "_log2_", libNameChIP, "_", libNameControl,
           "_featureFreq_chrPlot_winSize", winName, "_smooth", N,
           "_CSxRenan_step1Mb_IWGSCanalysis_MAPQ23_MAPQ2_MAPQ0_XM6_v", date, ".pdf"),
    height = 4, width = 10*length(chrName))
par(mfrow = c(1, length(chrName)))
par(mar = c(5.0, 9.0, 2.1, 9.0))
for(x in which(chrs %in% chrName)) {
  chrPartitionPlotCov3l_feature2(chrx = which(chrs %in% chrName),
                                 title = sub("c", "C", chrs[x]),
                                 cenStart = centromereStart[x],
                                 cenEnd = centromereEnd[x],
#                                 rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
#                                 rug1Col = "grey40",
                                 xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                 dat1A = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                 col1A = makeTransparent("darkorange"),
                                 dat1B = filt_chrProfilesChIPB[[x]]$filt_CPM,
                                 col1B = makeTransparent("darkgreen"),
                                 dat1C = filt_chrProfilesChIPC[[x]]$filt_CPM,
                                 col1C = makeTransparent("limegreen"),
                                 Ylab1 = bquote("Log"[2]*"("*.(sub("_\\w+", "", markChIP))*" ChIP/input)"),
                                 min1 = -max((minCPM*-1), maxCPM),
                                 max1 = maxCPM,
                                 legendLoc = "bottomright",
                                 legendLabs = c(expression("MAPQ ">="23"), expression("MAPQ ">="2"), expression("MAPQ ">="0"), "cM/Mb", "Genes"),
                                 xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                 dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                 col2A = makeTransparent("magenta2"),
                                 dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                 col2B = makeTransparent("dodgerblue"),
                                 Ylab2 = "",
                                 min2A = 0-maxFeatureA,
                                 max2A = maxFeatureA,
                                 min2B = 0-maxFeatureB,
                                 max2B = maxFeatureB)
}
dev.off()

pdf(paste0(plotDir, "Wheat",
           "_log2_", libNameChIP, "_", libNameControl,
           "_featureFreq_chrPlot_winSize", winName, "_smooth", N,
           "_CSxRenan_step1Mb_IWGSCanalysis_MAPQ23_MAPQ2_MAPQ0_XM6_v", date, ".pdf"),
    height = 4*7, width = 8*3)
par(mfrow = c(7, 3))
par(mar = c(5.0, 6.0, 2.1, 6.0))
for(x in which(chrs %in% chrs)) {
  chrPartitionPlotCov3l_feature2(chrx = which(chrs %in% chrs),
                                 title = sub("c", "C", chrs[x]),
                                 cenStart = centromereStart[x],
                                 cenEnd = centromereEnd[x],
#                                 rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
#                                 rug1Col = "grey40",
                                 xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                 dat1A = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                 col1A = makeTransparent("darkorange"),
                                 dat1B = filt_chrProfilesChIPB[[x]]$filt_CPM,
                                 col1B = makeTransparent("darkgreen"),
                                 dat1C = filt_chrProfilesChIPC[[x]]$filt_CPM,
                                 col1C = makeTransparent("limegreen"),
                                 Ylab1 = bquote("Log"[2]*"("*.(sub("_\\w+", "", markChIP))*" ChIP/input)"),
                                 min1 = -max((minCPM_chrs*-1), maxCPM_chrs),
                                 max1 = maxCPM_chrs,
                                 legendLoc = "bottomright",
                                 legendLabs = c(expression("MAPQ ">="23"), expression("MAPQ ">="2"), expression("MAPQ ">="0"), "cM/Mb", "Genes"),
                                 xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                 dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                 col2A = makeTransparent("magenta2"),
                                 dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                 col2B = makeTransparent("dodgerblue"),
                                 Ylab2 = "",
                                 min2A = 0-maxFeatureA_chrs,
                                 max2A = maxFeatureA_chrs,
                                 min2B = 0-maxFeatureB_chrs,
                                 max2B = maxFeatureB_chrs)
}
dev.off()

