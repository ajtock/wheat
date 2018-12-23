#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_ChIPseq_MNaseSeq.R H3K4me3 H3K4me3_Rep1_ChIP both 1Mb 1000000 15 forestgreen 281118 48

markChIP <- "H3K4me3"
libNameChIP <- "H3K4me3_Rep1_ChIP"
align <- "both"
winName <- "1Mb"
winSize <- 1000000
N <- 15
colour <- "forestgreen"
date <- 281118
cores <- 48

args <- commandArgs(trailingOnly = T)
markChIP <- args[1]
libNameChIP <- args[2]
align <- args[3]
winName <- args[4]
winSize <- as.numeric(args[5])
N <- as.numeric(args[6])
colour <- args[7]
date <- args[8]
cores <- as.numeric(args[9])

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)
library(plyr)
library(data.table)
library(varhandle)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
#chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
#chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,3])

plotDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/plots/"

## ChIP profile
covDirChIP <- paste0("/home/ajt200/analysis/wheat/",
                     markChIP, "/snakemake_ChIPseq/mapped/",
                     align, "/bg/")

profileChIP <- read.table(paste0(covDirChIP, libNameChIP, "_MappedOn_wheat_v1.0_lowXM_",
                                 align, "_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileChIP_bigWins <- profileChIP[profileChIP$V3-profileChIP$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileChIP <- profileChIP[profileChIP$V3-profileChIP$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end 
profileChIP_bigWinsList <- mclapply(seq_along(1:dim(profileChIP_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileChIP_bigWins[x,]$V2,
                      to = profileChIP_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileChIP_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIP_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileChIP_bigWins[x,]$V3)),
               V4 = as.numeric(profileChIP_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileChIP_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIP_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileChIP_bigWins[x,]$V4))
  }
}, mc.cores = cores)

profileChIP_bigWinsDT <- rbindlist(profileChIP_bigWinsList)
profileChIP <- rbind.fill(profileChIP, profileChIP_bigWinsDT)
profileChIP <- profileChIP[order(profileChIP$V1, profileChIP$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIP <- profileChIP[profileChIP$V1 == chrs[x],]
  if(chrProfileChIP[dim(chrProfileChIP)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIP[dim(chrProfileChIP)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIP[dim(chrProfileChIP)[1],]$V4))
  }
}, mc.cores = cores)
profileChIP_chrLenValsDT <- rbindlist(chrLenValsList)
profileChIP <- rbind.fill(profileChIP, profileChIP_chrLenValsDT)
profileChIP <- profileChIP[order(profileChIP$V1, profileChIP$V2),]

profileChIP <- data.frame(chr = as.character(profileChIP$V1),
                          window = as.integer(profileChIP$V2+1),
                          CPM = as.numeric(profileChIP$V4),
                          stringsAsFactors = F)

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  profileChIP[profileChIP$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfiles <- mclapply(seq_along(chrProfiles), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfiles[[x]]$CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfiles[[x]]$chr),
             window = as.integer(chrProfiles[[x]]$window),
             filt_CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfiles))

minCPM <- min(unlist(mclapply(seq_along(filt_chrProfiles),
  function(x) {
    min(filt_chrProfiles[[x]]$filt_CPM)
}, mc.cores = length(filt_chrProfiles))))
maxCPM <- max(unlist(mclapply(seq_along(filt_chrProfiles),
  function(x) {
    max(filt_chrProfiles[[x]]$filt_CPM)
}, mc.cores = length(filt_chrProfiles))))

# Plot
pdf(paste0(plotDir, "Wheat_", libNameChIP, "_",
           align, "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 22, width = 30)
par(mfrow = c(8, 3))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(x in 1:length(filt_chrProfiles)) {
  chrPlotCov(xplot = filt_chrProfiles[[x]]$window,
             title = chrs[x],
             cenStart = centromereStart[x],
             cenEnd = centromereEnd[x],
             dat1 = filt_chrProfiles[[x]]$filt_CPM,
             col1 = colour,
             Ylab1 = markChIP,
             min1 = minCPM, max1 = maxCPM)
}
dev.off()





## Define function to row-bind elements in bigWinsList
#rbind.fill.parallel <- function(bigWinsList) {
#  while(length(bigWinsList) > 1) {
#    bigWinsListIdx <- seq(from = 1, to = length(bigWinsList), by = 2)
#    bigWinsList <- mclapply(bigWinsListIdx, function(x) {
#      if(x==length(bigWinsList)) {
#        return(bigWinsList[[x]])
#      } else {
#        return(rbind.fill(bigWinsList[[x]], bigWinsList[[x+1]]))
#      }
#    }, mc.cores = 48)
#  }
#}
#
## Apply rbind.fill.parallel()
#rbind.fill.parallel(bigWinsList = profileChIP_bigWinsList)
