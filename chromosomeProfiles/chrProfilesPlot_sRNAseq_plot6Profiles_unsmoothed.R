#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_sRNAseq_plot6Profiles_unsmoothed.R CS+_2_LIB18613_LDI16228 both 1Mb 1000000 230118

libName1 <- "CS+_2_LIB18613_LDI16228"
align <- "both"
winName <- "1Mb"
winSize <- 1000000
date <- 150118

args <- commandArgs(trailingOnly = T)
libName1 <- args[1]
align <- args[2]
winName <- args[3]
winSize <- as.numeric(args[4])
date <- args[5]

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

names <- c("_20nt", "_21nt", "_22nt",
           "_23nt", "_24nt", "_33nt", "_34nt")
namesPlot <- c("20-nt", "21-nt", "22-nt",
               "23-nt", "24-nt", "33-nt", "34-nt")
colours <- c("red", "blue", "green2",
             "darkorange2", "purple3", "darkgreen", "deeppink")
makeTransparent <- function(thisColour, alpha = 150)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
coloursTransparent <- sapply(seq_along(colours), function(x) {
  makeTransparent(colours[x])
})

profilesDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/profiles/"
plotDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/plots/"

## Load total read counts for each sRNA size class
reads_per_sRNAsize <- sapply(seq_along(names), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/mapped/",
                    align, "/", libName1, "_MappedOn_wheat_v1.0_",
                    align, names[x],"_sort_reads.txt"))[,1]
})

# List of profiles, with one list element per sRNA size class
profiles <- lapply(seq_along(names), function(x) {
  read.table(paste0(profilesDir, "Wheat_sRNAseq_", libName1, "_",
                    align, names[x], "_chrPlot_winSize", winName, ".txt"),
             header = T)
})
profiles <- lapply(seq_along(profiles), function(x) {
  data.frame(chr = profiles[[x]]$chr,
             window = profiles[[x]]$window,
             CPM = profiles[[x]]$CPM*reads_per_sRNAsize[x])
})

# Separate into chromosomes
chrProfiles <- mclapply(seq_along(profiles), function(x) {
  mclapply(seq_along(chrs), function(y) {
    profiles[[x]][profiles[[x]] == chrs[y],]
  }, mc.cores = length(chrs))
}, mc.cores = length(profiles))


minCPM <- min(unlist(
  mclapply(seq_along(1:6), function(x) {
    min(unlist(
      mclapply(seq_along(chrs), function(y) {
        min(chrProfiles[[x]][[y]]$CPM)
      }, mc.cores = length(chrs))
    ))
}, mc.cores = 6)
))
maxCPM <- max(unlist(
  mclapply(seq_along(1:6), function(x) {
    max(unlist(
      mclapply(seq_along(chrs), function(y) {
        max(chrProfiles[[x]][[y]]$CPM)
      }, mc.cores = length(chrs))
    ))
}, mc.cores = 6)
))

minCPM7 <- min(unlist(
      mclapply(seq_along(chrs), function(y) {
        min(chrProfiles[[7]][[y]]$CPM)
      }, mc.cores = length(chrs))
))
maxCPM7 <- max(unlist(
      mclapply(seq_along(chrs), function(y) {
        max(chrProfiles[[7]][[y]]$CPM)
      }, mc.cores = length(chrs))
))

# Plot
pdf(paste0(plotDir, "Wheat_sRNAseq_", libName1, "_",
           align,"_7sRNAsizeClasses_chrPlot_winSize", winName,
           "_unsmoothed_v", date, ".pdf"),
    height = 22, width = 30)
par(mfrow = c(8, 3))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(y in 1:length(chrProfiles[[1]])) {
  chrPlotCov6(xplot = chrProfiles[[1]][[y]]$window,
              title = chrs[y],
              cenStart = centromereStart[y],
              cenEnd = centromereEnd[y],
              cols1 = coloursTransparent,
              Ylab1 = paste0("sRNAs"),
              min1 = minCPM, max1 = maxCPM,
              dat1A = chrProfiles[[1]][[y]]$CPM,
              dat1B = chrProfiles[[2]][[y]]$CPM,
              dat1C = chrProfiles[[3]][[y]]$CPM,
              dat1D = chrProfiles[[4]][[y]]$CPM,
              dat1E = chrProfiles[[5]][[y]]$CPM,
              dat1F = chrProfiles[[6]][[y]]$CPM,
              legendLoc = "topright",
              legendLabs = namesPlot[1:6])
}
plot.new()
plot.new()
for(y in 1:length(chrProfiles[[7]])) {
  chrPlotCov6(xplot = chrProfiles[[7]][[y]]$window,
              title = chrs[y],
              cenStart = centromereStart[y],
              cenEnd = centromereEnd[y],
              cols1 = coloursTransparent,
              Ylab1 = paste0("sRNAs"),
              min1 = minCPM7, max1 = maxCPM7,
              dat1A = chrProfiles[[7]][[y]]$CPM,
              dat1B = NULL,
              dat1C = NULL,
              dat1D = NULL,
              dat1E = NULL,
              dat1F = NULL,
              legendLoc = "topright",
              legendLabs = namesPlot[7])
}
dev.off()
