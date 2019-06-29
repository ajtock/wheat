#!/applications/R/R-3.5.0/bin/Rscript

# Plot library-size-normalized average coverage profiles around features and
# around equivalent random loci

# Usage:
# ./TEsuperfamProfilesPlot_ChIPseq_MNaseSeq_log2ChIPMnase.R both TEs 20bp 20 2kb 2000 141218

align <- "both"
featureName <- "TEs"
binName <- "20bp"
binSize <- 20
flankName <- "2kb"
flankSize <- 2000
date <- 141218

args <- commandArgs(trailingOnly = T)
align <- args[1]
featureName <- args[2]
binName <- args[3]
binSize <- as.numeric(args[4])
flankName <- args[5]
flankSize <- as.numeric(args[6])
date <- as.character(args[7])

# TE superfamily names
superfamName <- c(
                  "CACTA_DTC",
                  "Copia_LTR_RLC",
                  "Gypsy_LTR_RLG",
                  "Harbinger_DTH",
                  "hAT_DTA",
                  "Helitrons_DHH",
                  "LINE_RIX",
                  "Mariner_DTT",
                  "Mutator_DTM",
                  "SINE_SIX",
                  "Unclassified_class_2_DXX",
                  "Unclassified_LTR_RLX",
                  "Unclassified_repeats_XXX",
                  "Unclassified_with_TIRs_DTX"
                 )

plotDir <- sapply(seq_along(superfamName), function(x) {
  paste0("superfamily_", superfamName[x], "/")
})
sapply(seq_along(plotDir), function(x) {
  system(paste0("[ -d ", plotDir[x], " ] || mkdir ", plotDir[x]))
})

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)
library(doParallel)
registerDoParallel(cores = length(superfamName))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Sample names and directories
namesCombined <- c(
  "ASY1_Rep1_ChIP",
  "ASY1_Rep1_input",
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
  paste0("/home/ajt200/analysis/wheat/ASY1/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/ASY1/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/H2AZ/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/H3K27me1/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K36me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K9ac/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/CENH3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"),
  paste0("/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/TEProfiles/matrices/")
)

clrs <- c(
  "darkgreen",
  "seagreen",
  "dodgerblue", 
  "forestgreen",
  "green3",
  "magenta3",
  "firebrick1",
  "navy",
  "darkorange2",
  "green2",
  "deeppink",
  "grey60"
)

foreach(h = seq_along(superfamName)) %dopar% {
  
  # Load coverage matrices and mean of each column
  featureColMeansList <- mclapply(seq_along(namesCombined), function(x) {
    as.vector( colMeans(
                        read.table(paste0(inDirsCombined[x],
                                          namesCombined[x],
                                          "_MappedOn_wheat_v1.0_lowXM_",
                                          align, "_sort_norm_",
                                          superfamName[h], "_matrix_bin", binName,
                                          "_flank", flankName, ".tab"),
                                   skip = 3)
    , na.rm = T) )
  }, mc.cores = length(namesCombined))
  
  ranLocColMeansList <- mclapply(seq_along(namesCombined), function(x) {
    as.vector( colMeans(
                        read.table(paste0(inDirsCombined[x],
                                          namesCombined[x],
                                          "_MappedOn_wheat_v1.0_lowXM_",
                                          align, "_sort_norm_",
                                          superfamName[h], "_ranLoc_matrix_bin", binName,
                                          "_flank", flankName, ".tab"),
                                   skip = 3)
    , na.rm = T) )
  }, mc.cores = length(namesCombined))
  
  # Plot
  pdf(paste0(plotDir[h], "Wheat_", superfamName[h], "_", featureName, "_profiles_",
             align, "_bin", binName, "_flank", flankName,
             "_log2ChIPMnase_v", date, ".pdf"),
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

  print(h)

  for(x in 1:(length(namesCombined)-1)) {
    plotAvgCov(xplot = seq_along(featureColMeansList[[x]]),
               dat1 = featureColMeansList[[x]],
               ranDat1 = ranLocColMeansList[[x]],
               col1 = clrs[x],
               Ylab1 = namesCombined[x],
               flankSize = flankSize, binSize = binSize,
               flankLabL = "-2 kb", flankLabR = "+2 kb",
               featureStartLab = featureStartLab,
               featureEndLab = featureEndLab,
               ranLocStartLab = "Start",
               ranLocEndLab = "End")
    plotAvgCov(xplot = seq_along(featureColMeansList[[x]]),
               dat1 = log2((featureColMeansList[[x]]+1)/(featureColMeansList[[length(featureColMeansList)]]+1)),
               ranDat1 = log2((ranLocColMeansList[[x]]+1)/(ranLocColMeansList[[length(ranLocColMeansList)]]+1)),
               col1 = clrs[x],
               Ylab1 = bquote("Log"[2]*.(paste0("(", namesCombined[x], "/MNase)"))),
               flankSize = flankSize, binSize = binSize,
               flankLabL = "-2 kb", flankLabR = "+2 kb",
               featureStartLab = featureStartLab,
               featureEndLab = featureEndLab,
               ranLocStartLab = "Start",
               ranLocEndLab = "End")
  }
  plotAvgCov(xplot = seq_along(featureColMeansList[[length(featureColMeansList)]]),
             dat1 = featureColMeansList[[length(featureColMeansList)]],
             ranDat1 = ranLocColMeansList[[length(featureColMeansList)]],
             col1 = clrs[length(clrs)],
             Ylab1 = namesCombined[length(namesCombined)],
             flankSize = flankSize, binSize = binSize,
             flankLabL = "-2 kb", flankLabR = "+2 kb",
             featureStartLab = featureStartLab,
             featureEndLab = featureEndLab,
             ranLocStartLab = "Start",
             ranLocEndLab = "End")
  mtext(text = paste0(superfamName[h], " TEs"),
        outer = T, cex = 1, line = -1.5, at = 0.05)
  dev.off()
  
}
