#!/applications/R/R-4.0.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Change xblocks height to 5000 (ASY1) and 3000 (DMC1 and input) in chrPartitionPlotCov2_feature2 function

# Usage:
# ./chrProfilePlot_1homoeologuegroup_normNone_ChIP_MAPQ0_MAPQ2_MAPQ01_altColours.R ASY1_CS ASY1_CS_Rep1_ChIP ASY1_CS ASY1_CS_Rep1_ChIP ASY1_CS ASY1_CS_Rep1_ChIP both 1Mb 1000000 15 limegreen darkgreen blue 060621 'chr3A,chr3B,chr3D'
# ./chrProfilePlot_1homoeologuegroup_normNone_ChIP_MAPQ0_MAPQ2_MAPQ01_altColours.R DMC1 DMC1_Rep1_ChIP DMC1 DMC1_Rep1_ChIP DMC1 DMC1_Rep1_ChIP both 1Mb 1000000 15 limegreen darkgreen blue 060621 'chr3A,chr3B,chr3D'
# ./chrProfilePlot_1homoeologuegroup_normNone_ChIP_MAPQ0_MAPQ2_MAPQ01_altColours.R input H3_input_SRR6350669 input H3_input_SRR6350669 input H3_input_SRR6350669 both 1Mb 1000000 15 limegreen darkgreen blue 060621 'chr3A,chr3B,chr3D'

#markChIPA <- "DMC1"
#libNameChIPA <- "DMC1_Rep1_ChIP"
#markChIPB <- "DMC1"
#libNameChIPB <- "DMC1_Rep1_ChIP"
#markChIPC <- "DMC1"
#libNameChIPC <- "DMC1_Rep1_ChIP"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15
#colourA <- "limegreen"
#colourB <- "darkgreen"
#colourC <- "blue"
#date <- "060621"
#chrName <- unlist(strsplit("chr3A,chr3B,chr3D",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
markChIPB <- args[3]
libNameChIPB <- args[4]
markChIPC <- args[5]
libNameChIPC <- args[6]
align <- args[7]
winName <- args[8]
winSize <- as.numeric(args[9])
N <- as.numeric(args[10])
colourA <- args[11]
colourB <- args[12]
colourC <- args[13]
date <- args[14]
chrName <- unlist(strsplit(args[15],
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
colourC <- makeTransparent(colourC)

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
                       "CENH3_ChIP_SRR1686799",
                       "H3_input_SRR6350669")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPA, "/snakemake_ChIPseq_MAPQ0_XM6_normNone/mapped/",
                        align, "/bg/")
} else if(libNameChIPA %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                              "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                        markChIPA, "/snakemake_ChIPseq_MAPQ0_XM6_normNone/mapped/",
                        align, "/bg/")
} else {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPA, "/snakemake_ChIPseq_MAPQ0_XM6_normNone/mapped/",
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
                       "CENH3_ChIP_SRR1686799",
                       "H3_input_SRR6350669")) {
  covDirChIPB <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPB, "/snakemake_ChIPseq_MAPQ2_XM6_normNone/mapped/",
                        align, "/bg/")
} else if(libNameChIPB %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                              "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIPB <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                        markChIPB, "/snakemake_ChIPseq_MAPQ2_XM6_normNone/mapped/",
                        align, "/bg/")
} else {
  covDirChIPB <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPB, "/snakemake_ChIPseq_MAPQ2_XM6_normNone/mapped/",
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


## ChIPC profile
if(libNameChIPC %in% c("H3K4me3_ChIP_SRR6350668",
                       "H3K27me3_ChIP_SRR6350666",
                       "H3K36me3_ChIP_SRR6350670",
                       "H3K9ac_ChIP_SRR6350667",
                       "CENH3_ChIP_SRR1686799",
                       "H3_input_SRR6350669")) {
  covDirChIPC <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPC, "/snakemake_ChIPseq_MAPQ01_XM6_normNone/mapped/",
                        align, "/bg/")
} else if(libNameChIPC %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                              "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIPC <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                        markChIPC, "/snakemake_ChIPseq_MAPQ01_XM6_normNone/mapped/",
                        align, "/bg/")
} else {
  covDirChIPC <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPC, "/snakemake_ChIPseq_MAPQ01_XM6_normNone/mapped/",
                        align, "/bg/")
}
profileChIPC <- read.table(paste0(covDirChIPC, libNameChIPC, "_MappedOn_wheat_v1.0_lowXM_",
                                  align, "_sort_norm_binSize", winName, ".bedgraph"))
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

chrLenValsBList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPC <- profileChIPC[profileChIPC$V1 == chrs[x],]
  if(chrProfileChIPC[dim(chrProfileChIPC)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPC[dim(chrProfileChIPC)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPC[dim(chrProfileChIPC)[1],]$V4))
  }
}, mc.cores = detectCores())
profileChIPC_chrLenValsBDT <- rbindlist(chrLenValsBList)
profileChIPC <- rbind.fill(profileChIPC, profileChIPC_chrLenValsBDT)
profileChIPC <- profileChIPC[order(profileChIPC$V1, profileChIPC$V2),]

profileChIPC <- data.frame(chr = as.character(profileChIPC$V1),
                           window = as.integer(profileChIPC$V2+1),
                           CPM = as.numeric(profileChIPC$V4),
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

## ChIPC
# Get coverage within each window
#profileChIPC <- data.frame(chr = as.character(profileChIPC$chr),
#                           window = as.numeric(profileChIPC$window),
#                           CPM = as.numeric(profileChIPC$CPM),
#                           stringsAsFactors = F)
                        
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


minCPM_chrs <- min(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))
minCPM_chrs <- 0
maxCPM_chrs <- max(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))+1
}, mc.cores = length(filt_chrProfilesChIPA))))
#maxCPM_chrs <- 90

minCPM <- min(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    min(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesChIPA))))
minCPM <- 0
maxCPM <- max(unlist(mclapply(which(chrs %in% chrName),
  function(x) {
    max(c(filt_chrProfilesChIPA[[x]]$filt_CPM,
          filt_chrProfilesChIPB[[x]]$filt_CPM,
          filt_chrProfilesChIPC[[x]]$filt_CPM))+1
}, mc.cores = length(filt_chrProfilesChIPA))))
#maxCPM <- 90

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
           "_", libNameChIPA, "_MAPQ0_MAPQ2_MAPQ01",
           "_", align, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4, width = 8*length(chrName))
par(mfrow = c(1, length(chrName)))
#par(mar = c(5.0, 9.0, 2.1, 9.0))
par(mar = c(5.0, 6.0, 2.1, 6.0))
for(x in which(chrs %in% chrName)) {
  chrPartitionPlotCov3lCPM(chrx = which(chrs %in% chrName),
                           title = sub("c", "C", chrs[x]),
                           cenStart = centromereStart[x],
                           cenEnd = centromereEnd[x],
#                           rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
#                           rug1Col = "grey40",
                           xplot1 = filt_chrProfilesChIPA[[x]]$window,
                           dat1A = filt_chrProfilesChIPA[[x]]$filt_CPM,
                           col1A = colourA,
                           dat1B = filt_chrProfilesChIPB[[x]]$filt_CPM,
                           col1B = colourB,
                           dat1C = filt_chrProfilesChIPC[[x]]$filt_CPM,
                           col1C = colourC,
                           Ylab1 = paste0(sub("_\\w+", "", markChIPA), " coverage"),
                           min1 = minCPM,
                           max1 = maxCPM,
                           legendLoc = "right",
                           legendLabs = c(expression("MAPQ">=0),
                                          expression("MAPQ">=2),
                                          expression("MAPQ"<2))
                          )
}
dev.off()

pdf(paste0(plotDir, "Wheat",
           "_", libNameChIPA, "_MAPQ0_MAPQ2_MAPQ01",
           "_", align, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 4*7, width = 8*3)
par(mfrow = c(7, 3))
par(mar = c(5.0, 6.0, 2.1, 6.0))
for(x in which(chrs %in% chrs)) {
  chrPartitionPlotCov3lCPM(chrx = which(chrs %in% chrs),
                           title = sub("c", "C", chrs[x]),
                           cenStart = centromereStart[x],
                           cenEnd = centromereEnd[x],
#                           rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
#                           rug1Col = "grey40",
                           xplot1 = filt_chrProfilesChIPA[[x]]$window,
                           dat1A = filt_chrProfilesChIPA[[x]]$filt_CPM,
                           col1A = colourA,
                           dat1B = filt_chrProfilesChIPB[[x]]$filt_CPM,
                           col1B = colourB,
                           dat1C = filt_chrProfilesChIPC[[x]]$filt_CPM,
                           col1C = colourC,
                           Ylab1 = paste0(sub("_\\w+", "", markChIPA), " coverage"),
                           min1 = minCPM_chrs,
                           max1 = maxCPM_chrs,
                           legendLoc = "right",
                           legendLabs = c(expression("MAPQ">=0),
                                          expression("MAPQ">=2),
                                          expression("MAPQ"<2))
                          )
}
dev.off()
