#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_ChIPseq_IWGSC_Science_unsmoothed.R CENH3 CENH3_ChIP_SRR1686799 both 1Mb 1000000 deeppink 281118 48

markChIP <- "CENH3"
libNameChIP <- "CENH3_ChIP_SRR1686799"
align <- "both"
winName <- "1Mb"
winSize <- 1000000
colour <- "deeppink"
date <- 281118
cores <- 48

args <- commandArgs(trailingOnly = T)
markChIP <- args[1]
libNameChIP <- args[2]
align <- args[3]
winName <- args[4]
winSize <- as.numeric(args[5])
colour <- args[6]
date <- args[7]
cores <- as.numeric(args[8])

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
covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
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

minCPM <- min(unlist(mclapply(seq_along(chrProfiles),
  function(x) {
    min(chrProfiles[[x]]$CPM)
}, mc.cores = length(chrProfiles))))
maxCPM <- max(unlist(mclapply(seq_along(chrProfiles),
  function(x) {
    max(chrProfiles[[x]]$CPM)
}, mc.cores = length(chrProfiles))))

# Plot
pdf(paste0(plotDir, "Wheat_", libNameChIP, "_",
           align, "_chrPlot_winSize", winName,
           "_unsmoothed_v", date, ".pdf"),
    height = 22, width = 30)
par(mfrow = c(8, 3))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(x in 1:length(chrProfiles)) {
  chrPlotCov(xplot = chrProfiles[[x]]$window,
             title = chrs[x],
             cenStart = centromereStart[x],
             cenEnd = centromereEnd[x],
             dat1 = chrProfiles[[x]]$CPM,
             col1 = colour,
             Ylab1 = markChIP,
             min1 = minCPM, max1 = maxCPM)
}
dev.off()
