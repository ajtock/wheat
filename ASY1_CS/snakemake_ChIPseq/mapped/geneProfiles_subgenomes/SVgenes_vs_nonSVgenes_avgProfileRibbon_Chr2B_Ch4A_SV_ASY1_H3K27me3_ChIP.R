#!/applications/R/R-4.0.0/bin/Rscript

# Calculate and plot meta-profiles (means and 95% CIs) for genes
# within and outside of structural variant (SV) regions

# Chr2B SV region: 36,650,000-327,500,000
# Chr2B nonSV arm equivalent region: 400,000,000-750,000,000
# Chr4A SV region: 328,200,000-648,000,000
# Chr4A nonSV arm equivalent region: 50,000,000-225,000,000

# Usage:
# /applications/R/R-4.0.0/bin/Rscript SVgenes_vs_nonSVgenes_avgProfileRibbon_Chr2B_Ch4A_SV_ASY1_H3K27me3_ChIP.R both 3500 2000 2kb '2 kb' 20 20bp '0.02,0.96' 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' 

#align <- "both"
#bodyLength <- 3500
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
#                                        split = ",")))

args <- commandArgs(trailingOnly = T)
align <- args[1]
bodyLength <- as.numeric(args[2])
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
legendPos <- as.numeric(unlist(strsplit(args[8],
                                        split = ",")))
featureName <- unlist(strsplit(args[9],
                               split = ","))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

outDir <- "SVgenes_vs_nonSVgenes/"
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))


# Define plot titles
Chr2B_SV_featureNamePlot <- "Chr2B SV genes"
Chr2B_nonSV_featureNamePlot <- "Chr2B nonSV genes"
Chr2B_SV_ranLocNamePlot <- "Chr2B SV ranLoc"
Chr2B_nonSV_ranLocNamePlot <- "Chr2B nonSV ranLoc"
Chr4A_SV_featureNamePlot <- "Chr4A SV genes"
Chr4A_nonSV_featureNamePlot <- "Chr4A nonSV genes"
Chr4A_SV_ranLocNamePlot <- "Chr4A SV ranLoc"
Chr4A_nonSV_ranLocNamePlot <- "Chr4A nonSV ranLoc"

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]

## Load table of features grouped into quantiles
#if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
#  featuresDF <- read.table(paste0(outDir, "WesternEurope/features_", quantiles, "quantiles",
#                                  "_by_", libName, "_of_",
#                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
#                                  paste0(substring(featureName, first = 10, last = 16),
#                                         collapse = "_"), "_",
#                                  substring(featureName[1][1], first = 18), "_WesternEurope.txt"),
#                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
#} else {
#  featuresDF <- read.table(paste0(outDir, "WesternEurope/features_", quantiles, "quantiles",
#                                  "_by_", sub("_\\w+", "", libName), "_in_",
#                                  region, "_of_",
#                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
#                                  paste0(substring(featureName, first = 10, last = 16),
#                                         collapse = "_"), "_",
#                                  substring(featureName[1][1], first = 18), "_WesternEurope.txt"),
#                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
#}

# Load features table (which was used for generating the coverage matrices)
features <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                    paste0(substring(featureName[x], first = 10, last = 16),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 18), ".gff3"),
             header = F)
})
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature data.frames
if(length(featureName) == 3) {
 features <- do.call(rbind, features)
} else {
 features <- features[[1]]
}
colnames(features) <- c("chr", "version", "type", "start", "end", "score", "strand", "NA", "featureID")
#stopifnot(identical(as.character(featuresDF$featureID),
#                    as.character(features$V9)))

# Load ranLocs table (which was used for generating the coverage matrices)
ranLocs <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                    paste0(substring(featureName[x], first = 10, last = 16),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 18), "_randomLoci.bed"),
             header = F)
})
# If ranLocs from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature data.frames
if(length(featureName) == 3) {
 ranLocs <- do.call(rbind, ranLocs)
} else {
 ranLocs <- ranLocs[[1]]
}
colnames(ranLocs) <- c("chr", "start", "end", "ranLocID", "NA", "strand")
ranLocs$start <- as.integer(ranLocs$start+1)

# Get row indices for each feature group

# Chr2B SV region: 36,650,000-327,500,000
# Chr2B nonSV arm equivalent region: 400,000,000-750,000,000
Chr2B_SV_arm_region_start <- 36650000
Chr2B_SV_arm_region_end <- 327500000
Chr2B_nonSV_arm_region_start <- 400000000
Chr2B_nonSV_arm_region_end <- 750000000
# Chr4A SV region: 328,200,000-648,000,000
# Chr4A nonSV arm equivalent region: 50,000,000-225,000,000
Chr4A_SV_arm_region_start <- 328200000
Chr4A_SV_arm_region_end <- 648000000
Chr4A_nonSV_arm_region_start <- 50000000
Chr4A_nonSV_arm_region_end <- 225000000

Chr2B_SV_feature <- features[which(features$chr == "chr2B" &
                                   features$start >= Chr2B_SV_arm_region_start &
                                   features$end <= Chr2B_SV_arm_region_end),]
Chr2B_nonSV_feature <- features[which(features$chr == "chr2B" &
                                      features$start >= Chr2B_nonSV_arm_region_start &
                                      features$end <= Chr2B_nonSV_arm_region_end),]
Chr2B_SV_feature_indices <- which(features$chr == "chr2B" &
                                  features$start >= Chr2B_SV_arm_region_start &
                                  features$end <= Chr2B_SV_arm_region_end)
Chr2B_nonSV_feature_indices <- which(features$chr == "chr2B" &
                                     features$start >= Chr2B_nonSV_arm_region_start &
                                     features$end <= Chr2B_nonSV_arm_region_end)

Chr4A_SV_feature <- features[which(features$chr == "chr4A" &
                                   features$start >= Chr4A_SV_arm_region_start &
                                   features$end <= Chr4A_SV_arm_region_end),]
Chr4A_nonSV_feature <- features[which(features$chr == "chr4A" &
                                      features$start >= Chr4A_nonSV_arm_region_start &
                                      features$end <= Chr4A_nonSV_arm_region_end),]
Chr4A_SV_feature_indices <- which(features$chr == "chr4A" &
                                  features$start >= Chr4A_SV_arm_region_start &
                                  features$end <= Chr4A_SV_arm_region_end)
Chr4A_nonSV_feature_indices <- which(features$chr == "chr4A" &
                                     features$start >= Chr4A_nonSV_arm_region_start &
                                     features$end <= Chr4A_nonSV_arm_region_end)

Chr2B_SV_ranLoc <- ranLocs[which(ranLocs$chr == "chr2B" &
                                 ranLocs$start >= Chr2B_SV_arm_region_start &
                                 ranLocs$end <= Chr2B_SV_arm_region_end),]
Chr2B_nonSV_ranLoc <- ranLocs[which(ranLocs$chr == "chr2B" &
                                    ranLocs$start >= Chr2B_nonSV_arm_region_start &
                                    ranLocs$end <= Chr2B_nonSV_arm_region_end),]
Chr2B_SV_ranLoc_indices <- which(ranLocs$chr == "chr2B" &
                                 ranLocs$start >= Chr2B_SV_arm_region_start &
                                 ranLocs$end <= Chr2B_SV_arm_region_end)
Chr2B_nonSV_ranLoc_indices <- which(ranLocs$chr == "chr2B" &
                                    ranLocs$start >= Chr2B_nonSV_arm_region_start &
                                    ranLocs$end <= Chr2B_nonSV_arm_region_end)

Chr4A_SV_ranLoc <- ranLocs[which(ranLocs$chr == "chr4A" &
                                 ranLocs$start >= Chr4A_SV_arm_region_start &
                                 ranLocs$end <= Chr4A_SV_arm_region_end),]
Chr4A_nonSV_ranLoc <- ranLocs[which(ranLocs$chr == "chr4A" &
                                    ranLocs$start >= Chr4A_nonSV_arm_region_start &
                                    ranLocs$end <= Chr4A_nonSV_arm_region_end),]
Chr4A_SV_ranLoc_indices <- which(ranLocs$chr == "chr4A" &
                                 ranLocs$start >= Chr4A_SV_arm_region_start &
                                 ranLocs$end <= Chr4A_SV_arm_region_end)
Chr4A_nonSV_ranLoc_indices <- which(ranLocs$chr == "chr4A" &
                                    ranLocs$start >= Chr4A_nonSV_arm_region_start &
                                    ranLocs$end <= Chr4A_nonSV_arm_region_end)

### Random features
## Define function to randomly select n rows from
## a data.frame
#selectRandomFeatures <- function(features, n) {
#  return(features[sample(x = dim(features)[1],
#                         size = n,
#                         replace = FALSE),])
#}
#
## Define seed so that random selections are reproducible
#set.seed(93750174)

## Divide features into random sets of equal number,
## with the same number of genes per chromosome as
## above-defined libName-defined feature quantiles
#randomPCIndices <- lapply(1:quantiles, function(k) {
#  randomPCIndicesk <- NULL
#  for(i in 1:length(chrs)) {
#    randomPCfeatureskChr <- selectRandomFeatures(features = featuresDF[featuresDF$seqnames == chrs[i],],
#                                                 n = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k) &
#                                                                    featuresDF$seqnames == chrs[i],])[1])
#    randomPCIndicesk <- c(randomPCIndicesk, as.integer(rownames(randomPCfeatureskChr)))
#  }
#  randomPCIndicesk
#})
## Confirm per-chromosome feature numbers are the same for quantiles and random groupings
#lapply(seq_along(1:quantiles), function(k) {
#  sapply(seq_along(chrs), function(x) {
#    if(!identical(dim(featuresDF[randomPCIndices[[k]],][featuresDF[randomPCIndices[[k]],]$seqnames == chrs[x],]),
#                  dim(featuresDF[quantileIndices[[k]],][featuresDF[quantileIndices[[k]],]$seqnames == chrs[x],])))    {
#      stop("Quantile features and random features do not consist of the same number of features per chromosome")
#    }
#  })
#})


# Load feature matrices for each chromatin dataset and calculate log2(ChIP/control)
ChIPNames <- c(
               "Chr2B_SV_ASY1_Rep1_ChIP",
               "Chr4A_SV_ASY1_Rep1_ChIP",
               "Chr2B_SV_H3K27me3_Rep1_ChIP",
               "Chr4A_SV_H3K27me3_Rep1_ChIP",
               "ASY1_CS_Rep1_ChIP",
               "H3K27me3_ChIP_SRR6350666"
#               "DMC1_Rep1_ChIP"
              )
ChIPNamesDir <- c(
                  "Chr2B_SV_ASY1",
                  "Chr4A_SV_ASY1",
                  "Chr2B_SV_H3K27me3",
                  "Chr4A_SV_H3K27me3",
                  "ASY1_CS",
                  "H3K27me3"
#                  "DMC1"
                 )
ChIPNamesPlot <- c(
                   "Chr2B SV ASY1 (27 GB)",
                   "Chr4A SV ASY1 (23 GB)",
                   "Chr2B SV H3K27me3 (25 GB)",
                   "Chr4A SV H3K27me3 (10 GB)",
                   "WT ASY1 (43 GB)",
                   "WT H3K27me3 (53 GB)"
#                   "WT DMC1 (24 GB)"
                  )
ChIPColours <- c(
                 "mediumspringgreen",
                 "limegreen",
                 "magenta1",
                 "magenta4",
                 "darkgreen",
                 "firebrick1"
#                 "deepskyblue"
                )
log2ChIPNames <- ChIPNames
log2ChIPNamesPlot <- ChIPNamesPlot
log2ChIPColours <- ChIPColours
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  if(ChIPNames[x] %in% c("H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "CENH3_ChIP_SRR1686799")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else if(ChIPNames[x] %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                                "H3K27ac_Rep1_ChIP_SRR8126621")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    paste0("/home/ajt200/analysis/wheat/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  }
})

controlNames <- c(
                  "H3_input_SRR6350669",
                  "MNase_Rep1"
                 )
controlNamesDir <- c(
                     "input",
                     "MNase"
                    )
controlNamesPlot <- c(
                      "Input (25 GB)",
                      "MNase (70 GB)"
                     )
controlColours <- c(
                    "grey40",
                    "darkcyan"
                   )
controlDirs <- sapply(seq_along(controlNames), function(x) {
  if(controlNames[x] == "H3_input_SRR6350669") {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           controlNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else if(controlNames[x] == "MNase_Rep1") {
    paste0("/home/ajt200/analysis/wheat/",
           controlNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    if(!(controlNames %in% c("H3_input_SRR6350669", "MNase_Rep1"))) {
      stop(paste0("controlNames[", x, "] is neither H3_input_SRR6350669 nor MNase_Rep1"))
    }
  }
})


## control
# feature
control_featureMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x],
                                controlNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
control_featureMats <- mclapply(seq_along(control_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, control_featureMats[[x]])
  } else {
    control_featureMats[[x]][[1]]
  }
}, mc.cores = length(control_featureMats))

# ranLoc
control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x],
                                controlNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
control_ranLocMats <- mclapply(seq_along(control_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, control_ranLocMats[[x]])
  } else {
    control_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_ranLocMats))

## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, ChIP_featureMats[[x]])
  } else {
    ChIP_featureMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_featureMats))

ChIP_featureMatsChr2B <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP"
                        )) {
    ChIP_featureMats[[x]]/2
  } else {
    ChIP_featureMats[[x]]
  }
}, mc.cores = length(ChIP_featureMats))

ChIP_featureMatsChr4A <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP"
                        )) {
    ChIP_featureMats[[x]]/2
  } else {
    ChIP_featureMats[[x]]
  }
}, mc.cores = length(ChIP_featureMats))

# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, ChIP_ranLocMats[[x]])
  } else {
    ChIP_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_ranLocMats))

ChIP_ranLocMatsChr2B <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP"
                        )) {
    ChIP_ranLocMats[[x]]/2
  } else {
    ChIP_ranLocMats[[x]]
  }
}, mc.cores = length(ChIP_ranLocMats))

ChIP_ranLocMatsChr4A <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP"
                        )) {
    ChIP_ranLocMats[[x]]/2
  } else {
    ChIP_ranLocMats[[x]]
  }
}, mc.cores = length(ChIP_ranLocMats))

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP",
                         "ASY1_CS_Rep1_ChIP",
                         "DMC1_Rep1_ChIP",
                         "H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3K4me1_Rep1_ChIP_SRR8126618",
                         "H3K27ac_Rep1_ChIP_SRR8126621"
                        )) {
    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[2]]+1))
  }
}, mc.cores = length(ChIP_featureMats))

log2ChIP_featureMatsChr2B <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP",
                         "ASY1_CS_Rep1_ChIP",
                         "DMC1_Rep1_ChIP",
                         "H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3K4me1_Rep1_ChIP_SRR8126618",
                         "H3K27ac_Rep1_ChIP_SRR8126621"
                        )) {
    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr2B_SV_ASY1_Rep1_ChIP",
                           "Chr2B_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[1]]*2)+1))
    } else {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[1]])+1))
    }
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr2B_SV_ASY1_Rep1_ChIP",
                           "Chr2B_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[2]]*2)+1))
    } else {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[2]])+1))
    }
  }
}, mc.cores = length(ChIP_featureMats))

log2ChIP_featureMatsChr4A <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP",
                         "ASY1_CS_Rep1_ChIP",
                         "DMC1_Rep1_ChIP",
                         "H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3K4me1_Rep1_ChIP_SRR8126618",
                         "H3K27ac_Rep1_ChIP_SRR8126621"
                        )) {
    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr4A_SV_ASY1_Rep1_ChIP",
                           "Chr4A_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[1]]*2)+1))
    } else {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[1]])+1))
    }
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr4A_SV_ASY1_Rep1_ChIP",
                           "Chr4A_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[2]]*2)+1))
    } else {
      log2((ChIP_featureMats[[x]]+1)/((control_featureMats[[2]])+1))
    }
  }
}, mc.cores = length(ChIP_featureMats))

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP",
                         "ASY1_CS_Rep1_ChIP",
                         "DMC1_Rep1_ChIP",
                         "H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3K4me1_Rep1_ChIP_SRR8126618",
                         "H3K27ac_Rep1_ChIP_SRR8126621"
                        )) {
    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[2]]+1))
  }
}, mc.cores = length(ChIP_ranLocMats))

log2ChIP_ranLocMatsChr2B <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP",
                         "ASY1_CS_Rep1_ChIP",
                         "DMC1_Rep1_ChIP",
                         "H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3K4me1_Rep1_ChIP_SRR8126618",
                         "H3K27ac_Rep1_ChIP_SRR8126621"
                        )) {
    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr2B_SV_ASY1_Rep1_ChIP",
                           "Chr2B_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[1]]*2)+1))
    } else {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[1]])+1))
    }
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr2B_SV_ASY1_Rep1_ChIP",
                           "Chr2B_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[2]]*2)+1))
    } else {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[2]])+1))
    }
  }
}, mc.cores = length(ChIP_ranLocMats))

log2ChIP_ranLocMatsChr4A <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "Chr2B_SV_ASY1_Rep1_ChIP",
                         "Chr4A_SV_ASY1_Rep1_ChIP",
                         "Chr2B_SV_H3K27me3_Rep1_ChIP",
                         "Chr4A_SV_H3K27me3_Rep1_ChIP",
                         "ASY1_CS_Rep1_ChIP",
                         "DMC1_Rep1_ChIP",
                         "H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3K4me1_Rep1_ChIP_SRR8126618",
                         "H3K27ac_Rep1_ChIP_SRR8126621"
                        )) {
    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr4A_SV_ASY1_Rep1_ChIP",
                           "Chr4A_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[1]]*2)+1))
    } else {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[1]])+1))
    }
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    if(ChIPNames[x] %in% c(
                           "Chr4A_SV_ASY1_Rep1_ChIP",
                           "Chr4A_SV_H3K27me3_Rep1_ChIP"
                          )) {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[2]]*2)+1))
    } else {
      log2((ChIP_ranLocMats[[x]]+1)/((control_ranLocMats[[2]])+1))
    }
  }
}, mc.cores = length(ChIP_ranLocMats))


## log2ChIP
# Add column names
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_featureMatsChr2B[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMatsChr2B[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                               paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_featureMatsChr4A[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMatsChr4A[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                               paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the elements in the nested list correspond to coverage matrices for the different
# subsets of features and random loci
log2ChIP_mats <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  list(
       # Chr2B_SV_feature_indices
       log2ChIP_featureMatsChr2B[[x]][Chr2B_SV_feature_indices,],
       # Chr2B_nonSV_feature_indices
       log2ChIP_featureMats[[x]][Chr2B_nonSV_feature_indices,],
       # Chr2B_SV_ranLoc_indices
       log2ChIP_ranLocMatsChr2B[[x]][Chr2B_SV_ranLoc_indices,],
       # Chr2B_nonSV_ranLoc_indices
       log2ChIP_ranLocMats[[x]][Chr2B_nonSV_ranLoc_indices,],
       # Chr4A_SV_feature_indices
       log2ChIP_featureMatsChr4A[[x]][Chr4A_SV_feature_indices,],
       # Chr4A_nonSV_feature_indices
       log2ChIP_featureMats[[x]][Chr4A_nonSV_feature_indices,],
       # Chr4A_SV_ranLoc_indices
       log2ChIP_ranLocMatsChr4A[[x]][Chr4A_SV_ranLoc_indices,],
       # Chr4A_nonSV_ranLoc_indices
       log2ChIP_ranLocMats[[x]][Chr4A_nonSV_ranLoc_indices,]
      )
}, mc.cores = length(log2ChIP_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(log2ChIP_mats[[x]][[y]]),
               t(log2ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(log2ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  })
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    tidyDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    summaryDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
    summaryDFfeature_list_log2ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]])[1])
    summaryDFfeature_list_log2ChIP[[x]][[y]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_log2ChIP into
# a list of single data.frames containing all meta-profiles for plotting
Chr2B_SV_featureTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[1]]
})
Chr2B_nonSV_featureTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[2]]
})
Chr2B_SV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[3]]
})
Chr2B_nonSV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[4]]
})
Chr4A_SV_featureTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[5]]
})
Chr4A_nonSV_featureTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[6]]
})
Chr4A_SV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[7]]
})
Chr4A_nonSV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[8]]
})
names(Chr2B_SV_featureTmp) <- log2ChIPNamesPlot
names(Chr2B_nonSV_featureTmp) <- log2ChIPNamesPlot
names(Chr2B_SV_ranLocTmp) <- log2ChIPNamesPlot
names(Chr2B_nonSV_ranLocTmp) <- log2ChIPNamesPlot
names(Chr4A_SV_featureTmp) <- log2ChIPNamesPlot
names(Chr4A_nonSV_featureTmp) <- log2ChIPNamesPlot
names(Chr4A_SV_ranLocTmp) <- log2ChIPNamesPlot
names(Chr4A_nonSV_ranLocTmp) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- list(
  bind_rows(Chr2B_SV_featureTmp, .id = "libName"),
  bind_rows(Chr2B_nonSV_featureTmp, .id = "libName"),
  bind_rows(Chr2B_SV_ranLocTmp, .id = "libName"),
  bind_rows(Chr2B_nonSV_ranLocTmp, .id = "libName"),
  bind_rows(Chr4A_SV_featureTmp, .id = "libName"),
  bind_rows(Chr4A_nonSV_featureTmp, .id = "libName"),
  bind_rows(Chr4A_SV_ranLocTmp, .id = "libName"),
  bind_rows(Chr4A_nonSV_ranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
                                                   levels = log2ChIPNamesPlot)
}

# Define y-axis limits
ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
                       summaryDFfeature_log2ChIP[[2]]$CI_lower,
                       summaryDFfeature_log2ChIP[[3]]$CI_lower,
                       summaryDFfeature_log2ChIP[[4]]$CI_lower,
                       summaryDFfeature_log2ChIP[[5]]$CI_lower,
                       summaryDFfeature_log2ChIP[[6]]$CI_lower,
                       summaryDFfeature_log2ChIP[[7]]$CI_lower,
                       summaryDFfeature_log2ChIP[[8]]$CI_lower),
                 na.rm = T)
ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
                       summaryDFfeature_log2ChIP[[2]]$CI_upper,
                       summaryDFfeature_log2ChIP[[3]]$CI_upper,
                       summaryDFfeature_log2ChIP[[4]]$CI_upper,
                       summaryDFfeature_log2ChIP[[5]]$CI_upper,
                       summaryDFfeature_log2ChIP[[6]]$CI_upper,
                       summaryDFfeature_log2ChIP[[7]]$CI_upper,
                       summaryDFfeature_log2ChIP[[8]]$CI_upper),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = log2ChIPColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## Chr2B_SV_feature
summaryDFfeature <- summaryDFfeature_log2ChIP[[1]]
ggObj1_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_SV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_nonSV_feature
summaryDFfeature <- summaryDFfeature_log2ChIP[[2]]
ggObj2_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_nonSV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_SV_ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[3]]
ggObj3_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_SV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_nonSV_ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[4]]
ggObj4_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
#annotation_custom(legendLabs[[7]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_nonSV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_SV_feature
summaryDFfeature <- summaryDFfeature_log2ChIP[[5]]
ggObj5_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_SV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_nonSV_feature
summaryDFfeature <- summaryDFfeature_log2ChIP[[6]]
ggObj6_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_nonSV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_SV_ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[7]]
ggObj7_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_SV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_nonSV_ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[8]]
ggObj8_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/control)")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
#annotation_custom(legendLabs[[7]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_nonSV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_log2ChIP,
                                              ggObj2_combined_log2ChIP,
                                              ggObj3_combined_log2ChIP,
                                              ggObj4_combined_log2ChIP,
                                              ggObj5_combined_log2ChIP,
                                              ggObj6_combined_log2ChIP,
                                              ggObj7_combined_log2ChIP,
                                              ggObj8_combined_log2ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7,
                                                       8
                                                      ))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_avgProfiles_around_Chr2B_Chr4A_SV_nonSV_genes_ranLoc_",
              align, "_v020621.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*length(ggObjGA_combined), limitsize = FALSE)


## ChIP
# Add column names
for(x in seq_along(ChIP_featureMats)) {
  colnames(ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_featureMatsChr2B[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                            paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                            paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_ranLocMatsChr2B[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_featureMatsChr4A[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                            paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                            paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_ranLocMatsChr4A[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the elements in the nested list correspond to coverage matrices for the different
# subsets of features and random loci
ChIP_mats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  list(
       # Chr2B_SV_feature_indices
       ChIP_featureMatsChr2B[[x]][Chr2B_SV_feature_indices,],
       # Chr2B_nonSV_feature_indices
       ChIP_featureMats[[x]][Chr2B_nonSV_feature_indices,],
       # Chr2B_SV_ranLoc_indices
       ChIP_ranLocMatsChr2B[[x]][Chr2B_SV_ranLoc_indices,],
       # Chr2B_nonSV_ranLoc_indices
       ChIP_ranLocMats[[x]][Chr2B_nonSV_ranLoc_indices,],
       # Chr4A_SV_feature_indices
       ChIP_featureMatsChr4A[[x]][Chr4A_SV_feature_indices,],
       # Chr4A_nonSV_feature_indices
       ChIP_featureMats[[x]][Chr4A_nonSV_feature_indices,],
       # Chr4A_SV_ranLoc_indices
       ChIP_ranLocMatsChr4A[[x]][Chr4A_SV_ranLoc_indices,],
       # Chr4A_nonSV_ranLoc_indices
       ChIP_ranLocMats[[x]][Chr4A_nonSV_ranLoc_indices,]
      )
}, mc.cores = length(ChIP_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_ChIP <- mclapply(seq_along(ChIP_mats), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(ChIP_mats[[x]][[y]]),
               t(ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_ChIP  <- mclapply(seq_along(wideDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  })
}, mc.cores = length(wideDFfeature_list_ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    tidyDFfeature_list_ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_ChIP[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_ChIP  <- mclapply(seq_along(tidyDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_ChIP))

for(x in seq_along(summaryDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    summaryDFfeature_list_ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_ChIP[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
    summaryDFfeature_list_ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_ChIP[[x]][[y]])[1])
    summaryDFfeature_list_ChIP[[x]][[y]]$sem <- summaryDFfeature_list_ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_ChIP into
# a list of single data.frames containing all meta-profiles for plotting
Chr2B_SV_featureTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[1]]
})
Chr2B_nonSV_featureTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[2]]
})
Chr2B_SV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
Chr2B_nonSV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
Chr4A_SV_featureTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[5]]
})
Chr4A_nonSV_featureTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[6]]
})
Chr4A_SV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[7]]
})
Chr4A_nonSV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[8]]
})
names(Chr2B_SV_featureTmp) <- ChIPNamesPlot
names(Chr2B_nonSV_featureTmp) <- ChIPNamesPlot
names(Chr2B_SV_ranLocTmp) <- ChIPNamesPlot
names(Chr2B_nonSV_ranLocTmp) <- ChIPNamesPlot
names(Chr4A_SV_featureTmp) <- ChIPNamesPlot
names(Chr4A_nonSV_featureTmp) <- ChIPNamesPlot
names(Chr4A_SV_ranLocTmp) <- ChIPNamesPlot
names(Chr4A_nonSV_ranLocTmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(Chr2B_SV_featureTmp, .id = "libName"),
  bind_rows(Chr2B_nonSV_featureTmp, .id = "libName"),
  bind_rows(Chr2B_SV_ranLocTmp, .id = "libName"),
  bind_rows(Chr2B_nonSV_ranLocTmp, .id = "libName"),
  bind_rows(Chr4A_SV_featureTmp, .id = "libName"),
  bind_rows(Chr4A_nonSV_featureTmp, .id = "libName"),
  bind_rows(Chr4A_SV_ranLocTmp, .id = "libName"),
  bind_rows(Chr4A_nonSV_ranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                                   levels = ChIPNamesPlot)
}

# Define y-axis limits
ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
                       summaryDFfeature_ChIP[[2]]$CI_lower,
                       summaryDFfeature_ChIP[[3]]$CI_lower,
                       summaryDFfeature_ChIP[[4]]$CI_lower,
                       summaryDFfeature_ChIP[[5]]$CI_lower,
                       summaryDFfeature_ChIP[[6]]$CI_lower,
                       summaryDFfeature_ChIP[[7]]$CI_lower,
                       summaryDFfeature_ChIP[[8]]$CI_lower),
                 na.rm = T)
ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
                       summaryDFfeature_ChIP[[2]]$CI_upper,
                       summaryDFfeature_ChIP[[3]]$CI_upper,
                       summaryDFfeature_ChIP[[4]]$CI_upper,
                       summaryDFfeature_ChIP[[5]]$CI_upper,
                       summaryDFfeature_ChIP[[6]]$CI_upper,
                       summaryDFfeature_ChIP[[7]]$CI_upper,
                       summaryDFfeature_ChIP[[8]]$CI_upper),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = ChIPColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## Chr2B_SV_feature
summaryDFfeature <- summaryDFfeature_ChIP[[1]]
ggObj1_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_SV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_nonSV_feature
summaryDFfeature <- summaryDFfeature_ChIP[[2]]
ggObj2_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_nonSV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_SV_ranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[3]]
ggObj3_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_SV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_nonSV_ranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[4]]
ggObj4_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
#annotation_custom(legendLabs[[7]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_nonSV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_SV_feature
summaryDFfeature <- summaryDFfeature_ChIP[[5]]
ggObj5_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_SV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_nonSV_feature
summaryDFfeature <- summaryDFfeature_ChIP[[6]]
ggObj6_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_nonSV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_SV_ranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[7]]
ggObj7_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_SV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_nonSV_ranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[8]]
ggObj8_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "ChIP") +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
#annotation_custom(legendLabs[[7]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_nonSV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_ChIP,
                                              ggObj2_combined_ChIP,
                                              ggObj3_combined_ChIP,
                                              ggObj4_combined_ChIP,
                                              ggObj5_combined_ChIP,
                                              ggObj6_combined_ChIP,
                                              ggObj7_combined_ChIP,
                                              ggObj8_combined_ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7,
                                                       8
                                                      ))
ggsave(paste0(plotDir,
              "ChIP_avgProfiles_around_Chr2B_Chr4A_SV_nonSV_genes_ranLoc_",
              align, "_v020621.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*length(ggObjGA_combined), limitsize = FALSE)


## control
# Add column names
for(x in seq_along(control_featureMats)) {
  colnames(control_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the elements in the nested list correspond to coverage matrices for the different
# subsets of features and random loci
control_mats <- mclapply(seq_along(control_featureMats), function(x) {
  list(
       # Chr2B_SV_feature_indices
       control_featureMats[[x]][Chr2B_SV_feature_indices,],
       # Chr2B_nonSV_feature_indices
       control_featureMats[[x]][Chr2B_nonSV_feature_indices,],
       # Chr2B_SV_ranLoc_indices
       control_ranLocMats[[x]][Chr2B_SV_ranLoc_indices,],
       # Chr2B_nonSV_ranLoc_indices
       control_ranLocMats[[x]][Chr2B_nonSV_ranLoc_indices,],
       # Chr4A_SV_feature_indices
       control_featureMats[[x]][Chr4A_SV_feature_indices,],
       # Chr4A_nonSV_feature_indices
       control_featureMats[[x]][Chr4A_nonSV_feature_indices,],
       # Chr4A_SV_ranLoc_indices
       control_ranLocMats[[x]][Chr4A_SV_ranLoc_indices,],
       # Chr4A_nonSV_ranLoc_indices
       control_ranLocMats[[x]][Chr4A_nonSV_ranLoc_indices,]
      )
}, mc.cores = length(control_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_control <- mclapply(seq_along(control_mats), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = colnames(control_mats[[x]][[y]]),
               t(control_mats[[x]][[y]]))
  })
}, mc.cores = length(control_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_control[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  })
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    tidyDFfeature_list_control[[x]][[y]]$window <- factor(tidyDFfeature_list_control[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_control  <- mclapply(seq_along(tidyDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    summaryDFfeature_list_control[[x]][[y]]$window <- factor(summaryDFfeature_list_control[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
    summaryDFfeature_list_control[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]])[1])
    summaryDFfeature_list_control[[x]][[y]]$sem <- summaryDFfeature_list_control[[x]][[y]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]]$n-1)
    summaryDFfeature_list_control[[x]][[y]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
    summaryDFfeature_list_control[[x]][[y]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_control into
# a list of single data.frames containing all meta-profiles for plotting
Chr2B_SV_featureTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[1]]
})
Chr2B_nonSV_featureTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[2]]
})
Chr2B_SV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[3]]
})
Chr2B_nonSV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[4]]
})
Chr4A_SV_featureTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[5]]
})
Chr4A_nonSV_featureTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[6]]
})
Chr4A_SV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[7]]
})
Chr4A_nonSV_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[8]]
})
names(Chr2B_SV_featureTmp) <- controlNamesPlot
names(Chr2B_nonSV_featureTmp) <- controlNamesPlot
names(Chr2B_SV_ranLocTmp) <- controlNamesPlot
names(Chr2B_nonSV_ranLocTmp) <- controlNamesPlot
names(Chr4A_SV_featureTmp) <- controlNamesPlot
names(Chr4A_nonSV_featureTmp) <- controlNamesPlot
names(Chr4A_SV_ranLocTmp) <- controlNamesPlot
names(Chr4A_nonSV_ranLocTmp) <- controlNamesPlot
summaryDFfeature_control <- list(
  bind_rows(Chr2B_SV_featureTmp, .id = "libName"),
  bind_rows(Chr2B_nonSV_featureTmp, .id = "libName"),
  bind_rows(Chr2B_SV_ranLocTmp, .id = "libName"),
  bind_rows(Chr2B_nonSV_ranLocTmp, .id = "libName"),
  bind_rows(Chr4A_SV_featureTmp, .id = "libName"),
  bind_rows(Chr4A_nonSV_featureTmp, .id = "libName"),
  bind_rows(Chr4A_SV_ranLocTmp, .id = "libName"),
  bind_rows(Chr4A_nonSV_ranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_control)) {
  summaryDFfeature_control[[x]]$libName <- factor(summaryDFfeature_control[[x]]$libName,
                                                   levels = controlNamesPlot)
}

# Define y-axis limits
ymin_control <- min(c(summaryDFfeature_control[[1]]$CI_lower,
                       summaryDFfeature_control[[2]]$CI_lower,
                       summaryDFfeature_control[[3]]$CI_lower,
                       summaryDFfeature_control[[4]]$CI_lower,
                       summaryDFfeature_control[[5]]$CI_lower,
                       summaryDFfeature_control[[6]]$CI_lower,
                       summaryDFfeature_control[[7]]$CI_lower,
                       summaryDFfeature_control[[8]]$CI_lower),
                 na.rm = T)
ymax_control <- max(c(summaryDFfeature_control[[1]]$CI_upper,
                       summaryDFfeature_control[[2]]$CI_upper,
                       summaryDFfeature_control[[3]]$CI_upper,
                       summaryDFfeature_control[[4]]$CI_upper,
                       summaryDFfeature_control[[5]]$CI_upper,
                       summaryDFfeature_control[[6]]$CI_upper,
                       summaryDFfeature_control[[7]]$CI_upper,
                       summaryDFfeature_control[[8]]$CI_upper),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(controlNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(controlNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = controlColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## Chr2B_SV_feature
summaryDFfeature <- summaryDFfeature_control[[1]]
ggObj1_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[1]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[1]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[1]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_SV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_nonSV_feature
summaryDFfeature <- summaryDFfeature_control[[2]]
ggObj2_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[2]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[2]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[2]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_nonSV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_SV_ranLoc
summaryDFfeature <- summaryDFfeature_control[[3]]
ggObj3_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[3]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[3]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[3]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_SV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr2B_nonSV_ranLoc
summaryDFfeature <- summaryDFfeature_control[[4]]
ggObj4_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[4]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[4]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[4]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr2B_nonSV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_SV_feature
summaryDFfeature <- summaryDFfeature_control[[5]]
ggObj5_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[5]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[5]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[5]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_SV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_nonSV_feature
summaryDFfeature <- summaryDFfeature_control[[6]]
ggObj6_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[6]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[6]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "TSS",
                            "TTS",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[6]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_nonSV_featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_SV_ranLoc
summaryDFfeature <- summaryDFfeature_control[[7]]
ggObj7_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[7]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[7]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[7]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_SV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Chr4A_nonSV_ranLoc
summaryDFfeature <- summaryDFfeature_control[[8]]
ggObj8_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[8]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[8]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankNamePlot),
                            "Start",
                            "End",
                            paste0("+", flankNamePlot))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[8]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = "Control") +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(Chr4A_nonSV_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_control,
                                              ggObj2_combined_control,
                                              ggObj3_combined_control,
                                              ggObj4_combined_control,
                                              ggObj5_combined_control,
                                              ggObj6_combined_control,
                                              ggObj7_combined_control,
                                              ggObj8_combined_control
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7,
                                                       8
                                                      ))
ggsave(paste0(plotDir,
              "control_avgProfiles_around_Chr2B_Chr4A_SV_nonSV_genes_ranLoc_",
              align, "_v020621.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*length(ggObjGA_combined), limitsize = FALSE)
