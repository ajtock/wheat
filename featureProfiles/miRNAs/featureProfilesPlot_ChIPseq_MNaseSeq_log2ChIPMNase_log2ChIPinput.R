#!/applications/R/R-3.5.0/bin/Rscript

# Plot library-size-normalized average coverage profiles around features and
# around equivalent random loci

# Usage:
# ./featureProfilesPlot_ChIPseq_MNaseSeq_log2ChIPMNase_log2ChIPinput.R both miRNAs 5bp 5 500bp "500 bp" 500 220119

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
name1  <- "H3K9me2_Rep1_ChIP"
name2  <- "H3K4me3_Rep1_ChIP"
name3  <- "MNase_Rep1"
name4  <- "H3K4me3_ChIP_SRR6350668"
name5  <- "H3K27me3_ChIP_SRR6350666"
name6  <- "H3K36me3_ChIP_SRR6350670"
name7  <- "H3K9ac_ChIP_SRR6350667"
name8  <- "CENH3_ChIP_SRR1686799"
name9  <- "H3_input_SRR6350669"

inDir1 <- paste0("/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir2 <- paste0("/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir3 <- paste0("/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir4 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K4me3/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir5 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir6 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K36me3/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir7 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K9ac/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir8 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/CENH3/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")
inDir9 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/input/snakemake_ChIPseq/mapped/miRNAProfiles/matrices/")

namesCombined <-  c(
                    name1,
                    name2,
                    name3,
                    name4,
                    name5,
                    name6,
                    name7,
                    name8,
                    name9
                   )
inDirsCombined <- c(
                    inDir1,
                    inDir2,
                    inDir3,
                    inDir4,
                    inDir5,
                    inDir6,
                    inDir7,
                    inDir8,
                    inDir9
                   )

clrs <- c(
          "magenta3",
          "forestgreen",
          "grey40",
          "green3",
          "deepskyblue3",
          "darkorange",
          "dodgerblue3",
          "deeppink",
          "grey50"
         )

# Load coverage matrices and mean of each column
featureColMeanslist <- mclapply(seq_along(namesCombined), function(x) {
  as.vector( colMeans(
                      read.table(paste0(inDirsCombined[x],
                                        namesCombined[x],
                                        "_MappedOn_wheat_v1.0_lowXM_",
                                        align, "_sort_norm_",
                                        featureName, "_matrix_bin", binName,
                                        "_flank", flankName, ".tab"),
                                 skip = 3)
  , na.rm = T) )
}, mc.cores = length(namesCombined))

ranLocColMeanslist <- mclapply(seq_along(namesCombined), function(x) {
  as.vector( colMeans(
                      read.table(paste0(inDirsCombined[x],
                                        namesCombined[x],
                                        "_MappedOn_wheat_v1.0_lowXM_",
                                        align, "_sort_norm_",
                                        featureName, "_ranLoc_matrix_bin", binName,
                                        "_flank", flankName, ".tab"),
                                 skip = 3)
  , na.rm = T) )
}, mc.cores = length(namesCombined))

# Plot
pdf(paste0(plotDir, "Wheat_", featureName, "_profiles_",
           align, "_bin", binName, "_flank", flankName,
           "_log2ChIPinput_v", date, ".pdf"),
   height = 2.5*(length(namesCombined)), width = 3*4)
par(mfrow = c(length(namesCombined), 4))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

for(x in 1:2) {
  plotAvgCov(xplot = seq_along(featureColMeanslist[[x]]),
             dat1 = featureColMeanslist[[x]],
             ranDat1 = ranLocColMeanslist[[x]],
             col1 = clrs[x],
             Ylab1 = namesCombined[x],
             flankSize = flankSize, binSize = binSize,
             flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
  plotAvgCov(xplot = seq_along(featureColMeanslist[[x]]),
             dat1 = log2((featureColMeanslist[[x]]+1)/(featureColMeanslist[[3]]+1)),
             ranDat1 = log2((ranLocColMeanslist[[x]]+1)/(ranLocColMeanslist[[3]]+1)),
             col1 = clrs[x],
             Ylab1 = bquote("Log"[2]*.(paste0("(", namesCombined[x], "/MNase)"))),
             flankSize = flankSize, binSize = binSize,
             flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
}
for(x in 3:(length(featureColMeanslist)-1)) {
  plotAvgCov(xplot = seq_along(featureColMeanslist[[x]]),
             dat1 = featureColMeanslist[[x]],
             ranDat1 = ranLocColMeanslist[[x]],
             col1 = clrs[x],
             Ylab1 = namesCombined[x],
             flankSize = flankSize, binSize = binSize,
             flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
  plotAvgCov(xplot = seq_along(featureColMeanslist[[x]]),
             dat1 = log2((featureColMeanslist[[x]]+1)/(featureColMeanslist[[length(featureColMeanslist)]]+1)),
             ranDat1 = log2((ranLocColMeanslist[[x]]+1)/(ranLocColMeanslist[[length(featureColMeanslist)]]+1)),
             col1 = clrs[x],
             Ylab1 = bquote("Log"[2]*.(paste0("(", namesCombined[x], "/input)"))),
             flankSize = flankSize, binSize = binSize,
             flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
}
plotAvgCov(xplot = seq_along(featureColMeanslist[[8]]),
           dat1 = featureColMeanslist[[8]],
           ranDat1 = ranLocColMeanslist[[8]],
           col1 = clrs[8],
           Ylab1 = namesCombined[8],
           flankSize = flankSize, binSize = binSize,
           flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
           featureStartLab = featureStartLab,
           featureEndLab = featureEndLab,
           ranLocStartLab = "Start",
           ranLocEndLab = "End")
plotAvgCov(xplot = seq_along(featureColMeanslist[[8]]),
           dat1 = log2((featureColMeanslist[[8]]+1)/(featureColMeanslist[[3]]+1)),
           ranDat1 = log2((ranLocColMeanslist[[8]]+1)/(ranLocColMeanslist[[3]]+1)),
           col1 = clrs[8],
           Ylab1 = bquote("Log"[2]*.(paste0("(", namesCombined[8], "/MNase)"))),
           flankSize = flankSize, binSize = binSize,
           flankLabL = paste0("-", flankNamePlot), flankLabR = paste0("+", flankNamePlot),
           featureStartLab = featureStartLab,
           featureEndLab = featureEndLab,
           ranLocStartLab = "Start",
           ranLocEndLab = "End")
mtext(text = "miRNAs",
      outer = T, cex = 1, line = -1.5)
dev.off()
