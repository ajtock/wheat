#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_sRNAseq_unsmoothed.R CS+_2_LIB18613_LDI16228 both "_20nt" "20-nt" 1Mb 1000000 "grey30" 150118

libName1 <- "CS+_2_LIB18613_LDI16228"
align <- "both"
name1 <- "_20nt"
nameY1 <- "20-nt"
winName <- "1Mb"
winSize <- 1000000
colour1 <- "grey30"
date <- 150118

args <- commandArgs(trailingOnly = T)
libName1 <- args[1]
align <- args[2]
name1 <- args[3]
nameY1 <- args[4]
winName <- args[5]
winSize <- as.numeric(args[6])
colour1 <- args[7]
date <- args[8]

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

## profile1
covDir1 <- paste0("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/mapped/",
                  align, "/bg/")

profile1 <- read.table(paste0(covDir1, libName1, "_MappedOn_wheat_v1.0_",
                              align, name1, "_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profile1_bigWins <- profile1[profile1$V3-profile1$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profile1 <- profile1[profile1$V3-profile1$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end 
profile1_bigWinsList <- mclapply(seq_along(1:dim(profile1_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profile1_bigWins[x,]$V2,
                      to = profile1_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profile1_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profile1_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profile1_bigWins[x,]$V3)),
               V4 = as.numeric(profile1_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profile1_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profile1_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profile1_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profile1_bigWinsDT <- rbindlist(profile1_bigWinsList)
profile1 <- rbind.fill(profile1, profile1_bigWinsDT)
profile1 <- profile1[order(profile1$V1, profile1$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfile1 <- profile1[profile1$V1 == chrs[x],]
  if(chrProfile1[dim(chrProfile1)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfile1[dim(chrProfile1)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfile1[dim(chrProfile1)[1],]$V4))
  }
}, mc.cores = detectCores())
profile1_chrLenValsDT <- rbindlist(chrLenValsList)
profile1 <- rbind.fill(profile1, profile1_chrLenValsDT)
profile1 <- profile1[order(profile1$V1, profile1$V2),]

profile1 <- data.frame(chr = as.character(profile1$V1),
                       window = as.integer(profile1$V2+1),
                       CPM = as.numeric(profile1$V4),
                       stringsAsFactors = F)

chrProfiles1 <- mclapply(seq_along(chrs), function(x) {
  profile1[profile1$chr == chrs[x],]
}, mc.cores = length(chrs))

minCPM1 <- min(unlist(mclapply(seq_along(chrProfiles1),
  function(x) {
    min(chrProfiles1[[x]]$CPM)
}, mc.cores = length(chrProfiles1))))
maxCPM1 <- max(unlist(mclapply(seq_along(chrProfiles1),
  function(x) {
    max(chrProfiles1[[x]]$CPM)
}, mc.cores = length(chrProfiles1))))

# Plot
pdf(paste0(plotDir, "Wheat_sRNAseq_", libName1, "_",
           align, name1, "_chrPlot_winSize", winName,
           "_unsmoothed_v", date, ".pdf"),
    height = 22, width = 30)
par(mfrow = c(8, 3))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(x in 1:length(chrProfiles1)) {
  chrPlotCov(xplot = chrProfiles1[[x]]$window,
             title = chrs[x],
             cenStart = centromereStart[x],
             cenEnd = centromereEnd[x],
             dat1 = chrProfiles1[[x]]$CPM,
             col1 = colour1,
             Ylab1 = paste0(nameY1, " sRNAs"),
             min1 = minCPM1, max1 = maxCPM1)
}
dev.off()
