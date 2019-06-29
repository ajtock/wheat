#!/applications/R/R-3.5.0/bin/Rscript

# Plot library-size-normalized average coverage profiles around features and
# around equivalent random loci

# Usage:
# ./featureProfilesPlot_ChIPseq_MNaseSeq_log2ChIPMNase.R both genes 20bp 20 2kb 2000 080419

#align <- "both"
#featureName <- "genes"
#binName <- "20bp"
#binSize <- 20
#flankName <- "2kb"
#flankSize <- 2000
#date <- "080419"

args <- commandArgs(trailingOnly = T)
align <- args[1]
featureName <- args[2]
binName <- args[3]
binSize <- as.numeric(args[4])
flankName <- args[5]
flankSize <- as.numeric(args[6])
date <- as.character(args[7])

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)

plotDir <- "/home/ajt200/analysis/wheat/featureProfiles/genes/plots/"

# Sample names and directories
namesCombined <- c(
  "H2AZ_Rep1_ChIP",
  "H3K4me3_Rep1_ChIP",
  "H3K4me3_ChIP_SRR6350668",
  "H3K9me2_Rep1_ChIP",
  "H3K27me1_Rep1_ChIP",
  "H3K27me3_ChIP_SRR6350666",
  "H3K36me3_ChIP_SRR6350670",
  "H3K9ac_ChIP_SRR6350667",
  "CENH3_ChIP_SRR1686799",
  "MNase_Rep1"
)

inDirsCombined <- c(
  paste0("/home/ajt200/analysis/wheat/H2AZ/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/H3K27me1/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K36me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K9ac/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/CENH3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/geneProfiles/matrices/")
)

clrs <- c(
  "dodgerblue", 
  "forestgreen",
  "green3",
  "magenta3",
  "firebrick1",
  "navy",
  "darkorange2",
  "green2",
  "deeppink",
  "darkcyan"
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
           "_log2ChIPMNase_v", date, ".pdf"),
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

for(x in 1:(length(namesCombined)-1)) {
  plotAvgCov(xplot = seq_along(featureColMeanslist[[x]]),
             dat1 = featureColMeanslist[[x]],
             ranDat1 = ranLocColMeanslist[[x]],
             col1 = clrs[x],
             Ylab1 = namesCombined[x],
             flankSize = flankSize, binSize = binSize,
             flankLabL = "-2 kb", flankLabR = "+2 kb",
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
  plotAvgCov(xplot = seq_along(featureColMeanslist[[x]]),
             dat1 = log2((featureColMeanslist[[x]]+1)/(featureColMeanslist[[length(featureColMeanslist)]]+1)),
             ranDat1 = log2((ranLocColMeanslist[[x]]+1)/(ranLocColMeanslist[[length(ranLocColMeanslist)]]+1)),
             col1 = clrs[x],
             Ylab1 = bquote("Log"[2]*.(paste0("(", namesCombined[x], "/MNase)"))),
             flankSize = flankSize, binSize = binSize,
             flankLabL = "-2 kb", flankLabR = "+2 kb",
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
}
plotAvgCov(xplot = seq_along(featureColMeanslist[[length(featureColMeanslist)]]),
           dat1 = featureColMeanslist[[length(featureColMeanslist)]],
           ranDat1 = ranLocColMeanslist[[length(featureColMeanslist)]],
           col1 = clrs[length(clrs)],
           Ylab1 = namesCombined[length(namesCombined)],
           flankSize = flankSize, binSize = binSize,
           flankLabL = "-2 kb", flankLabR = "+2 kb",
           featureStartLab = featureStartLab,
           featureEndLab = featureEndLab,
           ranLocStartLab = "Start",
           ranLocEndLab = "End")
mtext(text = "Genes",
      outer = T, cex = 1, line = -1.5, at = 0.05)
dev.off()
