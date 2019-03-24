#!/applications/R/R-3.5.0/bin/Rscript

# Plot library-size-normalized average coverage profiles around features and
# around equivalent random loci

# Usage:
# ./featureProfilesPlot_sRNAseq.R both miRNAs 5bp 5 500bp "500 bp" 500 220119

#align <- "both"
#featureName <- "miRNAs"
#binName <- "5bp"
#binSize <- 5
#flankName <- "500bp"
#flankNamePlot <- "500 bp"
#flankSize <- 500
#date <- 220119

args <- commandArgs(trailingOnly = T)
align <- args[1]
featureName <- args[2]
binName <- args[3]
binSize <- as.numeric(args[4])
flankName <- args[5]
flankNamePlot <- args[6]
flankSize <- as.numeric(args[7])
date <- as.character(args[8])

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)

plotDir <- "/home/ajt200/analysis/wheat/featureProfiles/miRNAs/plots/"

# Sample names and directories

inDir <- "/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/mapped/miRNAProfiles/matrices/"
inDir_allSizes <- "/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/mapped/miRNAProfiles_allSizes/matrices/"

sizes <- c("_20nt", "_21nt", "_22nt", "_23nt",
           "_24nt", "_33nt", "_34nt")
sizesPlot <- c("20-nt", "21-nt", "22-nt", "23-nt",
               "24-nt", "33-nt", "34-nt")
clrs <- c("red", "blue", "green2", "darkorange2",
          "purple3", "darkgreen", "deeppink")

# Load coverage matrices and mean of each column
featureColMeanslist <- mclapply(seq_along(sizes), function(x) {
  as.vector( colMeans(
                      read.table(paste0(inDir[x],
                                        "CS+_2_LIB18613_LDI16228",
                                        "_MappedOn_wheat_v1.0_lowXM_",
                                        align, sizes[x], "_sort_norm_",
                                        featureName, "_matrix_bin", binName,
                                        "_flank", flankName, ".tab"),
                                 skip = 3)
  , na.rm = T) )
}, mc.cores = length(sizes))

ranLocColMeanslist <- mclapply(seq_along(sizes), function(x) {
  as.vector( colMeans(
                      read.table(paste0(inDir[x],
                                        "CS+_2_LIB18613_LDI16228",
                                        "_MappedOn_wheat_v1.0_lowXM_",
                                        align, sizes[x], "_sort_norm_",
                                        featureName, "_ranLoc_matrix_bin", binName,
                                        "_flank", flankName, ".tab"),
                                 skip = 3)
  , na.rm = T) )
}, mc.cores = length(sizes))

# Load sRNAseq (all sRNA sizes) coverage matrices and mean of each column
featureColMeans_allSizes <- as.vector( colMeans(
  read.table(paste0(inDir_allSizes[x],
                    "CS+_2_LIB18613_LDI16228",
                    "_MappedOn_wheat_v1.0_lowXM_",
                    align, "_sort_norm_",
                    featureName, "_matrix_bin", binName,
                    "_flank", flankName, ".tab"),
             skip = 3)
, na.rm = T) )

ranLocColMeans_allSizes <- as.vector( colMeans(
  read.table(paste0(inDir_allSizes[x],
                    "CS+_2_LIB18613_LDI16228",
                    "_MappedOn_wheat_v1.0_lowXM_",
                    align, "_sort_norm_",
                    featureName, "_ranLoc_matrix_bin", binName,
                    "_flank", flankName, ".tab"),
             skip = 3)
, na.rm = T) )

# Plot
pdf(paste0(plotDir, "Wheat_", featureName, "_profiles_",
           align, "_bin", binName, "_flank", flankName,
           "_sRNAseq_v", date, ".pdf"),
   height = 2.5*(length(sizes)+1), width = 3*2)
par(mfrow = c(length(sizes), 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

for(x in 1:length(featureColMeanslist)) {
  plotAvgCov(xplot = seq_along(featureColMeanslist[[x]]),
             dat1 = featureColMeanslist[[x]],
             ranDat1 = ranLocColMeanslist[[x]],
             col1 = clrs[x],
             Ylab1 = paste0(sizes[x], " sRNAs"),
             flankSize = flankSize, binSize = binSize,
             flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
}
plotAvgCov(xplot = seq_along(featureColMeans_allSizes),
           dat1 = featureColMeans_allSizes,
           ranDat1 = ranLocColMeans_allSizes,
           col1 = "grey50",
           Ylab1 = "All sRNAs",
           flankSize = flankSize, binSize = binSize,
           flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
           featureStartLab = featureStartLab,
           featureEndLab = featureEndLab,
           ranLocStartLab = "Start",
           ranLocEndLab = "End")
mtext(text = "miRNAs",
      outer = T, cex = 1, line = -1.5)
dev.off()
