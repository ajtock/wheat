#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_TEfamsFreq.R 1Mb 1000000 15 cyan3 291118 48 ./subfamilies/

winName <- "1Mb"
winSize <- 1000000
N <- 15
colour <- "cyan3"
date <- 291118
cores <- 48
inDir <- "./superfamilies/"

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
N <- as.numeric(args[3])
colour <- args[4]
date <- args[5]
cores <- as.numeric(args[6])
inDir <- args[7]

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)
library(doParallel)
library(stringr)
#library(plyr)
#library(data.table)
#library(varhandle)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,3])

plotDir <- paste0(inDir, "plots/")

# Read in feature frequency profiles and extract TE family names from file names
profileFiles <- system(paste0("ls ", inDir, "TE_frequency_per_", winName, "*.txt"), intern = T)
inDirCharLen <- length(unlist(strsplit(inDir, "")))
profileFilesCharLen <- sapply(seq_along(profileFiles), function(x) {
  length(unlist(strsplit(profileFiles[x], "")))
})
famNames <- sapply(seq_along(profileFiles), function(x) {
  substr(profileFiles[x],
         start = inDirCharLen + 34,
         # Assumes each file name has a 3-letter extension
         stop = profileFilesCharLen[x]-4)
})
if(grepl("subfamily", profileFiles[1]) == TRUE) {
  famNamesAbbrev <- sapply(seq_along(profileFiles), function(x) {
    str_match(profileFiles[x], "subfamily_(.*?).txt")[,2]
  })
} else {
  famNamesAbbrev <- famNames
}

profileFeatureList <- mclapply(seq_along(profileFiles), function(x) {
  read.table(profileFiles[x], header = T)
}, mc.cores = 48)

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

registerDoParallel(cores = cores)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(h = 1:length(profileFeatureList)) %dopar% {
  chrProfilesFeature <- lapply(seq_along(chrs), function(x) {
    profileFeatureList[[h]][profileFeatureList[[h]]$chr == chrs[x],]
  })

  filt_chrProfilesFeature <- lapply(seq_along(chrProfilesFeature), function(x) {
    filt_chrProfileFeature <- stats::filter(x = chrProfilesFeature[[x]]$winTEs,
                                            filter = f,
                                            sides = 2)
    filt_chrProfileFeature[1:flank] <- filt_chrProfileFeature[flank+1]
    filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank+1):length(filt_chrProfileFeature)] <- filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank)]
    data.frame(chr = as.character(chrProfilesFeature[[x]]$chr),
               window = as.integer(chrProfilesFeature[[x]]$window),
               filt_feature = as.numeric(filt_chrProfileFeature),
               stringsAsFactors = F)
  })

  minFeature <- min(unlist(lapply(seq_along(filt_chrProfilesFeature),
    function(x) {
      min(filt_chrProfilesFeature[[x]]$filt_feature)
  })))
  maxFeature <- max(unlist(lapply(seq_along(filt_chrProfilesFeature),
    function(x) {
      max(filt_chrProfilesFeature[[x]]$filt_feature)
  })))
 
  # Plot
  pdf(paste0(plotDir, "Wheat_", famNames[h],  
             "_chrPlot_winSize", winName, "_smooth", N,
             "_v", date, ".pdf"),
      height = 22, width = 30)
  par(mfrow = c(8, 3))
  par(mar = c(2.1, 4.1, 2.1, 4.1))
  par(mgp = c(3, 1, 0))
  
  for(x in 1:length(filt_chrProfilesFeature)) {
    chrPlotFeatureFreq(xplot = filt_chrProfilesFeature[[x]]$window,
                       title = chrs[x],
                       cenStart = centromereStart[x],
                       cenEnd = centromereEnd[x],
                       dat1 = filt_chrProfilesFeature[[x]]$filt_feature,
                       col1 = colour,
                       Ylab1 = famNamesAbbrev[h],
                       min1 = minFeature, max1 = maxFeature)
  }
  dev.off()
}
