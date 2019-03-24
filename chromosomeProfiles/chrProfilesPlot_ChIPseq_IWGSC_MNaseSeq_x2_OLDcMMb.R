#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_ChIPseq_IWGSC_MNaseSeq_x2_cMMb.R H3K36me3 H3K36me3_ChIP_SRR6350670 MNase MNase_Rep1 both 1Mb 1000000 15 darkorange2 darkcyan 280119

#markChIPA <- "H3K27me3"
#libNameChIPA <- "H3K27me3_ChIP_SRR6350666"
#markChIPB <- "H3K9ac"
#libNameChIPB <- "H3K9ac_ChIP_SRR6350667"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15
#colourA <- "navy"
#colourB <- "dodgerblue"
#date <- 280119

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
markChIPB <- args[3]
libNameChIPB <- args[4]
align <- args[5]
winName <- args[6]
winSize <- as.numeric(args[7])
N <- as.numeric(args[8])
colourA <- args[9]
colourB <- args[10]
date <- args[11]

makeTransparent <- function(thisColour, alpha = 150)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}

colourB <- makeTransparent(colourB)

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
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,3])

plotDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/plots/"

## ChIPA profile
covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                      markChIPA, "/snakemake_ChIPseq/mapped/",
                      align, "/bg/")

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

chrProfilesA <- mclapply(seq_along(chrs), function(x) {
  profileChIPA[profileChIPA$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfilesA <- mclapply(seq_along(chrProfilesA), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfilesA[[x]]$CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfilesA[[x]]$chr),
             window = as.integer(chrProfilesA[[x]]$window),
             filt_CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesA))

## ChIPB profile
covDirChIPB <- paste0("/home/ajt200/analysis/wheat/",
                      markChIPB, "/snakemake_ChIPseq/mapped/",
                      align, "/bg/")

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

chrProfilesB <- mclapply(seq_along(chrs), function(x) {
  profileChIPB[profileChIPB$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfilesB <- mclapply(seq_along(chrProfilesB), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfilesB[[x]]$CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfilesB[[x]]$chr),
             window = as.integer(chrProfilesB[[x]]$window),
             filt_CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesB))

# Set min and max values
minCPM <- min(unlist(mclapply(seq_along(filt_chrProfilesA),
  function(x) {
    min(c(filt_chrProfilesA[[x]]$filt_CPM,
          filt_chrProfilesB[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesA))))
maxCPM <- max(unlist(mclapply(seq_along(filt_chrProfilesA),
  function(x) {
    max(c(filt_chrProfilesA[[x]]$filt_CPM,
          filt_chrProfilesB[[x]]$filt_CPM))
}, mc.cores = length(filt_chrProfilesA))))


# Feature frequency chromosome profiles
profileFeatureFiles <- system("ls /home/irh25/wheat_files/chr*_w1mb_cMMb_v1.txt",
                              intern = T)
chrProfilesFeature <- mclapply(seq_along(chrs), function(x) {
  read.table(profileFeatureFiles[x], header = T)
}, mc.cores = length(chrs))

filt_chrProfilesFeature <- mclapply(seq_along(chrProfilesFeature), function(x) {
  filt_chrProfileFeature <- stats::filter(x = chrProfilesFeature[[x]]$win.vals,
                                          filter = f,
                                          sides = 2)
  filt_chrProfileFeature[1:flank] <- filt_chrProfileFeature[flank+1]
  filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank+1):length(filt_chrProfileFeature)] <- filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank)]
  data.frame(chr = as.character(chrProfilesFeature[[x]]$seqnames),
             window = as.integer(chrProfilesFeature[[x]]$start),
             filt_feature = as.numeric(filt_chrProfileFeature),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesFeature))

minFeature <- min(unlist(mclapply(seq_along(filt_chrProfilesFeature),
  function(x) {
    min(filt_chrProfilesFeature[[x]]$filt_feature)
}, mc.cores = length(filt_chrProfilesFeature))))
maxFeature <- max(unlist(mclapply(seq_along(filt_chrProfilesFeature),
  function(x) {
    max(filt_chrProfilesFeature[[x]]$filt_feature)
}, mc.cores = length(filt_chrProfilesFeature))))

# Plot
pdf(paste0(plotDir, "Wheat_", libNameChIPA, "_", libNameChIPB, "_",
           align, "_OLDcMMb_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 21, width = 30)
par(mfrow = c(7, 3))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(x in 1:length(filt_chrProfilesA)) {
  chrPlotCov2_cMMb(xplot1 = filt_chrProfilesA[[x]]$window,
                   xplot2 = filt_chrProfilesFeature[[x]]$window,
                   title = chrs[x],
                   cenStart = centromereStart[x],
                   cenEnd = centromereEnd[x],
                   dat1A = filt_chrProfilesA[[x]]$filt_CPM,
                   col1A = colourA,
                   dat1B = filt_chrProfilesB[[x]]$filt_CPM,
                   col1B = colourB,
                   Ylab1 = "Normalized coverage",
                   min1 = minCPM, max1 = maxCPM,
                   legendLoc = "top",
                   legendLabs = c(markChIPA, markChIPB), 
                   dat2 = filt_chrProfilesFeature[[x]]$filt_feature,
                   col2 = "red",
                   Ylab2 = "cM/Mb",
                   min2 = minFeature, max2 = maxFeature)
}
dev.off()
