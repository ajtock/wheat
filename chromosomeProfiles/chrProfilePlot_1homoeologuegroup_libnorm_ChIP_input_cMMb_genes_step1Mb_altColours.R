#!/applications/R/R-4.0.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Change xblocks height to 46.0 in chrPartitionPlotCov2_feature2 function

# Usage:
# ./chrProfilePlot_1homoeologuegroup_libnorm_ChIP_input_cMMb_genes_step1Mb_altColours.R DMC1 DMC1_Rep1_ChIP ASY1_CS ASY1_CS_Rep1_ChIP input H3_input_SRR6350669 input H3_input_SRR6350669 both 1Mb 1000000 15 limegreen darkgreen 120421 '_MAPQ0_XM6' 'chr3A,chr3B,chr3D'

#markChIPA <- "DMC1"
#libNameChIPA <- "DMC1_Rep1_ChIP"
#markChIPB <- "ASY1_CS"
#libNameChIPB <- "ASY1_CS_Rep1_ChIP"
#markControlA <- "input"
#libNameControlA <- "H3_input_SRR6350669"
#markControlB <- "input"
#libNameControlB <- "H3_input_SRR6350669"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15
#colourA <- "limegreen"
#colourB <- "darkgreen"
#date <- "120421"
#params <- "_MAPQ0_XM6" 
#chrName <- unlist(strsplit("chr3A,chr3B,chr3D",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
markChIPB <- args[3]
libNameChIPB <- args[4]
markControlA <- args[5]
libNameControlA <- args[6]
markControlB <- args[7]
libNameControlB <- args[8]
align <- args[9]
winName <- args[10]
winSize <- as.numeric(args[11])
N <- as.numeric(args[12])
colourA <- args[13]
colourB <- args[14]
date <- args[15]
params <- args[16]
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
colourB <- makeTransparent(colourB)

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
                        markChIPA, "/snakemake_ChIPseq", params, "/mapped/",
                        align, "/bg/")
}
profileChIPA <- read.table(paste0(covDirChIPA, libNameChIPA, "_MappedOn_wheat_v1.0_lowXM_",
                                  align, "_sort_norm_binSize", winName, ".bedgraph"))
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

## ChIPB profile
if(libNameChIPB %in% c("H3K4me3_ChIP_SRR6350668",
                       "H3K27me3_ChIP_SRR6350666",
                       "H3K36me3_ChIP_SRR6350670",
                       "H3K9ac_ChIP_SRR6350667",
                       "CENH3_ChIP_SRR1686799")) {
  covDirChIPB <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPB, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else if(libNameChIPB %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                              "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIPB <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                        markChIPB, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else {
  covDirChIPB <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPB, "/snakemake_ChIPseq", params, "/mapped/",
                        align, "/bg/")
}
profileChIPB <- read.table(paste0(covDirChIPB, libNameChIPB, "_MappedOn_wheat_v1.0_lowXM_",
                                  align, "_sort_norm_binSize", winName, ".bedgraph"))
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

chrLenValsBList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPB <- profileChIPB[profileChIPB$V1 == chrs[x],]
  if(chrProfileChIPB[dim(chrProfileChIPB)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPB[dim(chrProfileChIPB)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPB[dim(chrProfileChIPB)[1],]$V4))
  }
}, mc.cores = detectCores())
profileChIPB_chrLenValsBDT <- rbindlist(chrLenValsBList)
profileChIPB <- rbind.fill(profileChIPB, profileChIPB_chrLenValsBDT)
profileChIPB <- profileChIPB[order(profileChIPB$V1, profileChIPB$V2),]

profileChIPB <- data.frame(chr = as.character(profileChIPB$V1),
                           window = as.integer(profileChIPB$V2+1),
                           CPM = as.numeric(profileChIPB$V4),
                           stringsAsFactors = F)

## ControlA profile
if(libNameControlA == "MNase_Rep1") {
  covDirControlA <- paste0("/home/ajt200/analysis/wheat/",
                           "MNase/snakemake_ChIPseq/mapped/", align, "/bg/")
  profileControlA <- read.table(paste0(covDirControlA, "MNase_Rep1_MappedOn_wheat_v1.0_lowXM_",
                                       align, "_sort_norm_binSize", winName, ".bedgraph"))
} else if(libNameControlA == "H3_input_SRR6350669") {
  covDirControlA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                           "input/snakemake_ChIPseq", params, "/mapped/", align, "/bg/")
  profileControlA <- read.table(paste0(covDirControlA, "H3_input_SRR6350669_MappedOn_wheat_v1.0_lowXM_",
                                       align, "_sort_norm_binSize", winName, ".bedgraph"))
} else {
  if(!(libNameControlA %in% c("MNase_Rep1", "H3_input_SRR6350669"))) {
    stop("libNameControlA is neither MNase_Rep1 nor H3_input_SRR6350669")
  }
}
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

## ControlB profile
if(libNameControlB == "MNase_Rep1") {
  covDirControlB <- paste0("/home/ajt200/analysis/wheat/",
                           "MNase/snakemake_ChIPseq/mapped/", align, "/bg/")
  profileControlB <- read.table(paste0(covDirControlB, "MNase_Rep1_MappedOn_wheat_v1.0_lowXM_",
                                       align, "_sort_norm_binSize", winName, ".bedgraph"))
} else if(libNameControlB == "H3_input_SRR6350669") {
  covDirControlB <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                           "input/snakemake_ChIPseq", params, "/mapped/", align, "/bg/")
  profileControlB <- read.table(paste0(covDirControlB, "H3_input_SRR6350669_MappedOn_wheat_v1.0_lowXM_",
                                       align, "_sort_norm_binSize", winName, ".bedgraph"))
} else {
  if(!(libNameControlB %in% c("MNase_Rep1", "H3_input_SRR6350669"))) {
    stop("libNameControlB is neither MNase_Rep1 nor H3_input_SRR6350669")
  }
}
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

chrLenValsBList <- mclapply(seq_along(chrs), function (x) {
  chrProfileControl <- profileControlB[profileControlB$V1 == chrs[x],]
  if(chrProfileControl[dim(chrProfileControl)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileControl[dim(chrProfileControl)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileControl[dim(chrProfileControl)[1],]$V4))
  }
}, mc.cores = detectCores())
profileControlB_chrLenValsBDT <- rbindlist(chrLenValsBList)
profileControlB <- rbind.fill(profileControlB, profileControlB_chrLenValsBDT)
profileControlB <- profileControlB[order(profileControlB$V1, profileControlB$V2),]

profileControlB <- data.frame(chr = as.character(profileControlB$V1),
                           window = as.integer(profileControlB$V2+1),
                           CPM = as.numeric(profileControlB$V4),
                           stringsAsFactors = F)

## ChIPA
# Get coverage within each window
#profileChIPA <- data.frame(chr = as.character(profileChIPA$chr),
#                           window = as.numeric(profileChIPA$window),
#                           CPM = as.numeric(profileChIPA$CPM),
#                           stringsAsFactors = F)
                        
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

## ChIPB
# Get coverage within each window
#profileChIPB <- data.frame(chr = as.character(profileChIPB$chr),
#                           window = as.numeric(profileChIPB$window),
#                           CPM = as.numeric(profileChIPB$CPM),
#                           stringsAsFactors = F)
                        
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

## ControlA
# Get coverage within each window
#profileControlA <- data.frame(chr = as.character(profileControlA$chr),
#                              window = as.numeric(profileControlA$window),
#                              CPM = as.numeric(profileControlA$CPM),
#                              stringsAsFactors = F)
                        
chrProfilesControlA <- mclapply(seq_along(chrs), function(x) {
  profileControlA[profileControlA$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfilesControlA <- mclapply(seq_along(chrProfilesControlA), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfilesControlA[[x]]$CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfilesControlA[[x]]$chr),
             window = as.integer(chrProfilesControlA[[x]]$window),
             filt_CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesControlA))


minCPM_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesControlA[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))
#minCPM_chrs <- 40
maxCPM_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesControlA[[x]]$filt_CPM))+1
}, mc.cores = length(filt_chrProfilesChIPA))))
#maxCPM_chrs <- 100

minCPM <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesControlA[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))
#minCPM <- 40
maxCPM <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesControlA[[x]]$filt_CPM))+1
}, mc.cores = length(filt_chrProfilesChIPA))))
#maxCPM <- 100

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
           "_", libNameChIPA, "_", libNameChIPB, "_", libNameControlA,
           "_", align, "_featureFreq_chrPlot_winSize", winName, "_smooth", N,
           "_CSxRenan_step1Mb_IWGSCanalysis_v", date, params, ".pdf"),
    height = 4, width = 8*length(chrName))
par(mfrow = c(1, length(chrName)))
#par(mar = c(5.0, 9.0, 2.1, 9.0))
par(mar = c(5.0, 6.0, 2.1, 6.0))
for(x in which(chrs %in% chrName)) {
  chrPartitionPlotCov3lCPM_feature2(chrx = which(chrs %in% chrName),
                                 title = sub("c", "C", chrs[x]),
                                 cenStart = centromereStart[x],
                                 cenEnd = centromereEnd[x],
#                                 rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
#                                 rug1Col = "grey40",
                                 xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                 dat1A = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                 col1A = colourA,
                                 dat1B = filt_chrProfilesChIPB[[x]]$filt_CPM,
                                 col1B = colourB,
                                 dat1C = filt_chrProfilesControlA[[x]]$filt_CPM,
                                 col1C = "red",
                                 Ylab1 = bquote(.(toupper(align)) ~ .(substring(params, 2)) ~ "(CPM)"),
                                 min1 = minCPM,
                                 max1 = maxCPM,
                                 legendLoc = "topright",
                                 legendLabs = c(sub("_\\w+", "", markChIPA), sub("_\\w+", "", markChIPB), "Input", "cM/Mb", "Genes"),
                                 xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                 dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                 col2A = "cyan2",
                                 dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                 col2B = "grey80",
                                 Ylab2 = "",
                                 min2A = minFeatureA,
                                 max2A = maxFeatureA,
                                 min2B = minFeatureB,
                                 max2B = maxFeatureB)
}
dev.off()

pdf(paste0(plotDir, "Wheat",
           "_", libNameChIPA, "_", libNameChIPB, "_", libNameControlA,
           "_", align, "_featureFreq_chrPlot_winSize", winName, "_smooth", N,
           "_CSxRenan_step1Mb_IWGSCanalysis_v", date, params, ".pdf"),
    height = 4*7, width = 8*3)
par(mfrow = c(7, 3))
par(mar = c(5.0, 6.0, 2.1, 6.0))
for(x in which(chrs %in% chrs)) {
  chrPartitionPlotCov3lCPM_feature2(chrx = which(chrs %in% chrs),
                                 title = sub("c", "C", chrs[x]),
                                 cenStart = centromereStart[x],
                                 cenEnd = centromereEnd[x],
#                                 rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
#                                 rug1Col = "grey40",
                                 xplot1 = filt_chrProfilesChIPA[[x]]$window,
                                 dat1A = filt_chrProfilesChIPA[[x]]$filt_CPM,
                                 col1A = colourA,
                                 dat1B = filt_chrProfilesChIPB[[x]]$filt_CPM,
                                 col1B = colourB,
                                 dat1C = filt_chrProfilesControlA[[x]]$filt_CPM,
                                 col1C = "red",
                                 Ylab1 = bquote(.(toupper(align)) ~ .(substring(params, 2)) ~ "(CPM)"),
                                 min1 = minCPM_chrs,
                                 max1 = maxCPM_chrs,
                                 legendLoc = "topright",
                                 legendLabs = c(sub("_\\w+", "", markChIPA), sub("_\\w+", "", markChIPB), "Input", "cM/Mb", "Genes"),
                                 xplot2 = filt_chrProfilesFeatureA[[x]]$window,
                                 dat2A = filt_chrProfilesFeatureA[[x]]$filt_feature,
                                 col2A = "cyan2",
                                 dat2B = filt_chrProfilesFeatureB[[x]]$filt_feature,
                                 col2B = "grey80",
                                 Ylab2 = "",
                                 min2A = 0,
                                 max2A = maxFeatureA_chrs,
                                 min2B = 0,
                                 max2B = maxFeatureB_chrs)
}
dev.off()

