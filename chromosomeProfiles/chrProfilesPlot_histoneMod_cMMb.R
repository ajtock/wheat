#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_histoneMod_cMMb.R MNase MNase_Rep1 both 1Mb 1000000 15 darkcyan 100619 200

#markChIPA <- "MNase"
#libNameChIPA <- "MNase_Rep1"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15
#colourA <- "darkcyan"
#date <- 100619
#minMarkerDist <- 200

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
align <- args[3]
winName <- args[4]
winSize <- as.numeric(args[5])
N <- as.numeric(args[6])
colourA <- args[7]
date <- args[8]
minMarkerDist <- as.numeric(args[9])

makeTransparent <- function(thisColour, alpha = 200)
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

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])

plotDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/plots/"

## ChIPA profile
if(libNameChIPA %in% c("H3K4me3_ChIP_SRR6350668",
                       "H3K27me3_ChIP_SRR6350666",
                       "H3K36me3_ChIP_SRR6350670",
                       "H3K9ac_ChIP_SRR6350667",
                       "CENH3_ChIP_SRR1686799")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
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

minCPM <- min(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    min(filt_chrProfilesChIPA[[x]]$filt_CPM)
}, mc.cores = length(filt_chrProfilesChIPA))))
maxCPM <- max(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    max(filt_chrProfilesChIPA[[x]]$filt_CPM)
}, mc.cores = length(filt_chrProfilesChIPA))))


# Feature frequency chromosome profiles
cMMbProfile <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                                 as.character(minMarkerDist), "bp_", winName, ".txt"),
                          header = T)
chrProfilesFeature <- mclapply(seq_along(chrs), function(x) {
  cMMbProfile[cMMbProfile$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfilesFeature <- mclapply(seq_along(chrProfilesFeature), function(x) {
  filt_chrProfileFeature <- stats::filter(x = chrProfilesFeature[[x]]$cMMb,
                                          filter = f,
                                          sides = 2)
  filt_chrProfileFeature[1:flank] <- filt_chrProfileFeature[flank+1]
  filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank+1):length(filt_chrProfileFeature)] <- filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank)]
  data.frame(chr = as.character(chrProfilesFeature[[x]]$chr),
             window = as.integer(chrProfilesFeature[[x]]$window),
             filt_feature = as.numeric(filt_chrProfileFeature),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesFeature))

minFeature <- min(unlist(mclapply(seq_along(filt_chrProfilesFeature),
  function(x) {
    min(filt_chrProfilesFeature[[x]]$filt_feature, na.rm = T)
}, mc.cores = length(filt_chrProfilesFeature))))
maxFeature <- max(unlist(mclapply(seq_along(filt_chrProfilesFeature),
  function(x) {
    max(filt_chrProfilesFeature[[x]]$filt_feature, na.rm = T)
}, mc.cores = length(filt_chrProfilesFeature))))


# Plot
pdf(paste0(plotDir, "Wheat_", libNameChIPA, "_",
           align, "_cMMb_chrPlot_winSize", winName, "_smooth", N,
           "_minInterMarkerDist", as.character(minMarkerDist),
           "bp_v", date, ".pdf"),
    height = 21, width = 30)
par(mfrow = c(7, 3))
par(mar = c(2.1, 4.5, 2.1, 4.5))
par(mgp = c(3, 1, 0))

for(x in 1:length(filt_chrProfilesChIPA)) {
  chrPlotCov_cMMb(xplot1 = filt_chrProfilesChIPA[[x]]$window,
                  xplot2 = filt_chrProfilesFeature[[x]]$window,
                  title = chrs[x],
                  cenStart = centromereStart[x],
                  cenEnd = centromereEnd[x],
                  dat1 = filt_chrProfilesChIPA[[x]]$filt_CPM,
                  col1 = colourA,
                  Ylab1 = "Normalized coverage",
                  min1 = 0,
                  max1 = maxCPM,
                  legendLoc = "top",
                  legendLabs = c(markChIPA),
                  dat2 = filt_chrProfilesFeature[[x]]$filt_feature,
                  col2 = "cyan",
                  Ylab2 = "cM/Mb",
                  min2 = 0, max2 = maxFeature)
}
dev.off()
