#!/applications/R/R-3.4.0/bin/Rscript

# Plot average coverage profiles with 95% CIs around
# Cluster 1 genes annotated with the  "defense response" GO term; i.e.,
# clusters_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters/GO/cluster1_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters_of_genes_in_Agenome_genomewide_GO_BP/cluster1_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters_of_genes_in_Agenome_genomewide_GO_BP_enrichment_GO\:0006952.txt

# Usage:
# /applications/R/R-3.4.0/bin/Rscript defense_response_genes_avgProfileRibbon.R genes_in_Agenome_genomewide 'Defense_response' 3500 2000 2kb '2 kb' 20 20bp promoters '1' '4' '0006952' both ASY1_CS_Rep1_ChIP ASY1_CS purple4

featureName <- "genes_in_Agenome_genomewide"
featureNamePlot <- "Defense_response"
bodyLength <- 3500
upstream <- 2000
downstream <- 2000
flankName <- "2kb"
flankNamePlot <- "2 kb"
binSize <- 20
binName <- "20bp"
region <- "promoters"
clusterNo <- "1"
clusterLast <- "4"
GO_ID <- "0006952"
align <- "both"

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
featureNamePlot <- args[2]
bodyLength <- as.numeric(args[3])
upstream <- as.numeric(args[4])
downstream <- as.numeric(args[5])
flankName <- args[6]
flankNamePlot <- args[7]
binSize <- as.numeric(args[8])
binName <- args[9]
region <- args[10]
clusterNo <- as.character(args[11])
clusterLast <- as.character(args[12])
GO_ID <- as.character(args[13])
align <- as.character(args[14])

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

outDir <- paste0("clusters_by_log2_", libName,
                 "_control_in_", region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

IDs <- as.character(read.table(paste0("clusters_by_log2_ASY1_CS_Rep1_ChIP_control_in_", region,
                                      "/GO/cluster", clusterNo, "_of_", clusterLast,
                                      "_by_log2_ASY1_CS_Rep1_ChIP_control_in_", region, "_of_", featureName,
                                      "_GO_BP/cluster", clusterNo, "_of_", clusterLast,
                                      "_by_log2_ASY1_CS_Rep1_ChIP_control_in_", region, "_of_", featureName,
                                      "_GO_BP_enrichment_GO:", GO_ID, ".txt"),
                               colClasses = c("NULL", NA), header = F)$V2)
IDs <- unlist(strsplit(x = IDs,
                       split = ","))

# Load features
features <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                              substring(featureName, first = 10), ".gff3"),
                       header = F)
featureIDs <- sub(pattern = "\\.\\d+", replacement = "",
                  features$V9)
ID_indices <- which(featureIDs %in% IDs)

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
# and sort by decreasing log2mat1RegionRowMeans
ChIPNames <- c(
               "ASY1_CS_Rep1_ChIP",
               "DMC1_Rep1_ChIP",
               "H2AZ_Rep1_ChIP",
               "H3K4me3_Rep1_ChIP",
               "H3K4me1_Rep1_ChIP_SRR8126618",
               "H3K27ac_Rep1_ChIP_SRR8126621",
               "H3K27me3_ChIP_SRR6350666",
               "H3K9me2_Rep1_ChIP",
               "H3K27me1_Rep1_ChIP"
              )
ChIPNamesDir <- c(
                  "ASY1_CS",
                  "DMC1",
                  "H2AZ",
                  "H3K4me3",
                  "H3K4me1",
                  "H3K27ac",
                  "H3K27me3",
                  "H3K9me2",
                  "H3K27me1"
                 )
ChIPNamesPlot <- c(
                   "ASY1",
                   "DMC1",
                   "H2A.Z",
                   "H3K4me3",
                   "H3K4me1",
                   "H3K27ac",
                   "H3K27me3",
                   "H3K9me2",
                   "H3K27me1"
                  )
ChIPColours <- c(
                 "purple4",
                 "green2",
                 "dodgerblue",
                 "forestgreen",
                 "goldenrod1",
                 "orange",
                 "navy",
                 "magenta3",
                 "firebrick1"
                )
otherNames <- c(
                "MNase_Rep1",
                "DNaseI_Rep1_SRR8447247"
               )
otherNamesDir <- c(
                   "MNase",
                   "DNaseI"
                  )
otherNamesPlot <- c(
                    "MNase",
                    "DNaseI"
                   )
otherColours <- c(
                  "darkcyan",
                  "purple"
                 )
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
otherDirs <- sapply(seq_along(otherNames), function(x) {
  if(otherNames[x] %in% c("MNase_Rep1")) {
    paste0("/home/ajt200/analysis/wheat/",
           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else if(otherNames[x] %in% c("DNaseI_Rep1_SRR8447247")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    if(!(otherNames %in% c("MNase_Rep1", "DNaseI_Rep1_SRR8447247"))) {
      stop("otherNames[x] is neither MNase_Rep1 nor DNaseI_Rep1_SRR8447247")
    }
  }
})
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x],
                              ChIPNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))
other_featureMats <- mclapply(seq_along(otherNames), function(x) {
  as.matrix(read.table(paste0(otherDirs[x],
                              otherNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(otherNames))
# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x],
                              ChIPNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_ranLoc_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))
other_ranLocMats <- mclapply(seq_along(otherNames), function(x) {
  as.matrix(read.table(paste0(otherDirs[x],
                              otherNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_ranLoc_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(otherNames))

# Load control matrices
controlNames <- c(
                  "H3_input_SRR6350669",
                  "MNase_Rep1"
                 )
controlNamesDir <- c(
                     "input",
                     "MNase"
                    )
controlNamesPlot <- c(
                      "Input",
                      "MNase"
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
      stop("controlNames[x] is neither H3_input_SRR6350669 nor MNase_Rep1")
    }
  }
})
# feature
control_featureMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))
# ranLoc
control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_ranLoc_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
# feature
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(ChIPNames[x] %in% c(
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

# ranLoc
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(ChIPNames[x] %in% c(
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

# Add column names, and
# extract only row numbers (features and ranLoc) in ID_indices
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}
log2ChIP_featureMats <- lapply(seq_along(log2ChIP_featureMats), function(x) {
  log2ChIP_featureMats[[x]][ID_indices,]
})
log2ChIP_ranLocMats <- lapply(seq_along(log2ChIP_ranLocMats), function(x) {
  log2ChIP_ranLocMats[[x]][ID_indices,]
})
for(x in seq_along(other_featureMats)) {
  colnames(other_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                        paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(other_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                       paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                       paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}
other_featureMats <- lapply(seq_along(other_featureMats), function(x) {
  other_featureMats[[x]][ID_indices,]
})
other_ranLocMats <- lapply(seq_along(other_ranLocMats), function(x) {
  other_ranLocMats[[x]][ID_indices,]
})
for(x in seq_along(control_featureMats)) {
  colnames(control_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                         paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                         paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}
control_featureMats <- lapply(seq_along(control_featureMats), function(x) {
  control_featureMats[[x]][ID_indices,]
})
control_ranLocMats <- lapply(seq_along(control_ranLocMats), function(x) {
  control_ranLocMats[[x]][ID_indices,]
})

## feature
# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  data.frame(window = colnames(log2ChIP_featureMats[[x]]),
             t(log2ChIP_featureMats[[x]]))
}, mc.cores = length(log2ChIP_featureMats))

wideDFfeature_list_other <- mclapply(seq_along(other_featureMats), function(x) {
  data.frame(window = colnames(other_featureMats[[x]]),
             t(other_featureMats[[x]]))
}, mc.cores = length(other_featureMats))

wideDFfeature_list_control <- mclapply(seq_along(control_featureMats), function(x) {
  data.frame(window = colnames(control_featureMats[[x]]),
             t(control_featureMats[[x]]))
}, mc.cores = length(control_featureMats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  gather(data  = wideDFfeature_list_log2ChIP[[x]],
         key   = feature,
         value = coverage,
         -window)
}, mc.cores = length(wideDFfeature_list_log2ChIP))

tidyDFfeature_list_other  <- mclapply(seq_along(wideDFfeature_list_other), function(x) {
  gather(data  = wideDFfeature_list_other[[x]],
         key   = feature,
         value = coverage,
         -window)
}, mc.cores = length(wideDFfeature_list_other))

tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  gather(data  = wideDFfeature_list_control[[x]],
         key   = feature,
         value = coverage,
         -window)
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  tidyDFfeature_list_log2ChIP[[x]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]]$window,
                                                    levels = as.character(wideDFfeature_list_log2ChIP[[x]]$window))
}

for(x in seq_along(tidyDFfeature_list_other)) {
  tidyDFfeature_list_other[[x]]$window <- factor(tidyDFfeature_list_other[[x]]$window,
                                                 levels = as.character(wideDFfeature_list_other[[x]]$window))
}

for(x in seq_along(tidyDFfeature_list_control)) {
  tidyDFfeature_list_control[[x]]$window <- factor(tidyDFfeature_list_control[[x]]$window,
                                                   levels = as.character(wideDFfeature_list_control[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]]$window),
             n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFfeature_list_log2ChIP[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFfeature_list_log2ChIP[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFfeature_list_log2ChIP[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  summaryDFfeature_list_log2ChIP[[x]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]]$window,
                                                       levels = as.character(wideDFfeature_list_log2ChIP[[x]]$window))
  summaryDFfeature_list_log2ChIP[[x]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]])[1])
  summaryDFfeature_list_log2ChIP[[x]]$sem <- summaryDFfeature_list_log2ChIP[[x]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]]$n-1)
  summaryDFfeature_list_log2ChIP[[x]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]]$mean -
    qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]]$n-1)*summaryDFfeature_list_log2ChIP[[x]]$sem
  summaryDFfeature_list_log2ChIP[[x]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]]$mean +
    qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]]$n-1)*summaryDFfeature_list_log2ChIP[[x]]$sem
}

names(summaryDFfeature_list_log2ChIP) <- ChIPNamesPlot

summaryDFfeature_list_other  <- mclapply(seq_along(tidyDFfeature_list_other), function(x) {
  data.frame(window = as.character(wideDFfeature_list_other[[x]]$window),
             n      = tapply(X     = tidyDFfeature_list_other[[x]]$coverage,
                             INDEX = tidyDFfeature_list_other[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFfeature_list_other[[x]]$coverage,
                             INDEX = tidyDFfeature_list_other[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFfeature_list_other[[x]]$coverage,
                             INDEX = tidyDFfeature_list_other[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFfeature_list_other))

for(x in seq_along(summaryDFfeature_list_other)) {
  summaryDFfeature_list_other[[x]]$window <- factor(summaryDFfeature_list_other[[x]]$window,
                                                    levels = as.character(wideDFfeature_list_other[[x]]$window))
  summaryDFfeature_list_other[[x]]$winNo <- factor(1:dim(summaryDFfeature_list_other[[x]])[1])
  summaryDFfeature_list_other[[x]]$sem <- summaryDFfeature_list_other[[x]]$sd/sqrt(summaryDFfeature_list_other[[x]]$n-1)
  summaryDFfeature_list_other[[x]]$CI_lower <- summaryDFfeature_list_other[[x]]$mean -
    qt(0.975, df = summaryDFfeature_list_other[[x]]$n-1)*summaryDFfeature_list_other[[x]]$sem
  summaryDFfeature_list_other[[x]]$CI_upper <- summaryDFfeature_list_other[[x]]$mean +
    qt(0.975, df = summaryDFfeature_list_other[[x]]$n-1)*summaryDFfeature_list_other[[x]]$sem
}

names(summaryDFfeature_list_other) <- otherNamesPlot

summaryDFfeature_list_control  <- mclapply(seq_along(tidyDFfeature_list_control), function(x) {
  data.frame(window = as.character(wideDFfeature_list_control[[x]]$window),
             n      = tapply(X     = tidyDFfeature_list_control[[x]]$coverage,
                             INDEX = tidyDFfeature_list_control[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFfeature_list_control[[x]]$coverage,
                             INDEX = tidyDFfeature_list_control[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFfeature_list_control[[x]]$coverage,
                             INDEX = tidyDFfeature_list_control[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  summaryDFfeature_list_control[[x]]$window <- factor(summaryDFfeature_list_control[[x]]$window,
                                                      levels = as.character(wideDFfeature_list_control[[x]]$window))
  summaryDFfeature_list_control[[x]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]])[1])
  summaryDFfeature_list_control[[x]]$sem <- summaryDFfeature_list_control[[x]]$sd/sqrt(summaryDFfeature_list_control[[x]]$n-1)
  summaryDFfeature_list_control[[x]]$CI_lower <- summaryDFfeature_list_control[[x]]$mean -
    qt(0.975, df = summaryDFfeature_list_control[[x]]$n-1)*summaryDFfeature_list_control[[x]]$sem
  summaryDFfeature_list_control[[x]]$CI_upper <- summaryDFfeature_list_control[[x]]$mean +
    qt(0.975, df = summaryDFfeature_list_control[[x]]$n-1)*summaryDFfeature_list_control[[x]]$sem
}

names(summaryDFfeature_list_control) <- controlNamesPlot

# Convert list summaryDFfeature_list_log2ChIP into a single data.frame for plotting
summaryDFfeature_log2ChIP <- bind_rows(summaryDFfeature_list_log2ChIP, .id = "libName")
summaryDFfeature_log2ChIP$libName <- factor(summaryDFfeature_log2ChIP$libName,
                                            levels = names(summaryDFfeature_list_log2ChIP))

summaryDFfeature_other <- bind_rows(summaryDFfeature_list_other, .id = "libName")
summaryDFfeature_other$libName <- factor(summaryDFfeature_other$libName,
                                         levels = names(summaryDFfeature_list_other))

summaryDFfeature_control <- bind_rows(summaryDFfeature_list_control, .id = "libName")
summaryDFfeature_control$libName <- factor(summaryDFfeature_control$libName,
                                           levels = names(summaryDFfeature_list_control))

## ranLoc
# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFranLoc_list_log2ChIP <- mclapply(seq_along(log2ChIP_ranLocMats), function(x) {
  data.frame(window = colnames(log2ChIP_ranLocMats[[x]]),
             t(log2ChIP_ranLocMats[[x]]))
}, mc.cores = length(log2ChIP_ranLocMats))

wideDFranLoc_list_other <- mclapply(seq_along(other_ranLocMats), function(x) {
  data.frame(window = colnames(other_ranLocMats[[x]]),
             t(other_ranLocMats[[x]]))
}, mc.cores = length(other_ranLocMats))

wideDFranLoc_list_control <- mclapply(seq_along(control_ranLocMats), function(x) {
  data.frame(window = colnames(control_ranLocMats[[x]]),
             t(control_ranLocMats[[x]]))
}, mc.cores = length(control_ranLocMats))

# Convert into tidy data.frame (long format)
tidyDFranLoc_list_log2ChIP  <- mclapply(seq_along(wideDFranLoc_list_log2ChIP), function(x) {
  gather(data  = wideDFranLoc_list_log2ChIP[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list_log2ChIP))

tidyDFranLoc_list_other  <- mclapply(seq_along(wideDFranLoc_list_other), function(x) {
  gather(data  = wideDFranLoc_list_other[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list_other))

tidyDFranLoc_list_control  <- mclapply(seq_along(wideDFranLoc_list_control), function(x) {
  gather(data  = wideDFranLoc_list_control[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFranLoc_list_log2ChIP)) {
  tidyDFranLoc_list_log2ChIP[[x]]$window <- factor(tidyDFranLoc_list_log2ChIP[[x]]$window,
                                                   levels = as.character(wideDFranLoc_list_log2ChIP[[x]]$window))
}

for(x in seq_along(tidyDFranLoc_list_other)) {
  tidyDFranLoc_list_other[[x]]$window <- factor(tidyDFranLoc_list_other[[x]]$window,
                                                levels = as.character(wideDFranLoc_list_other[[x]]$window))
}

for(x in seq_along(tidyDFranLoc_list_control)) {
  tidyDFranLoc_list_control[[x]]$window <- factor(tidyDFranLoc_list_control[[x]]$window,
                                                  levels = as.character(wideDFranLoc_list_control[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (ranLocs) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFranLoc_list_log2ChIP  <- mclapply(seq_along(tidyDFranLoc_list_log2ChIP), function(x) {
  data.frame(window = as.character(wideDFranLoc_list_log2ChIP[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_log2ChIP[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_log2ChIP[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranLoc_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_log2ChIP[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranLoc_list_log2ChIP))

for(x in seq_along(summaryDFranLoc_list_log2ChIP)) {
  summaryDFranLoc_list_log2ChIP[[x]]$window <- factor(summaryDFranLoc_list_log2ChIP[[x]]$window,
                                                      levels = as.character(wideDFranLoc_list_log2ChIP[[x]]$window))
  summaryDFranLoc_list_log2ChIP[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list_log2ChIP[[x]])[1])
  summaryDFranLoc_list_log2ChIP[[x]]$sem <- summaryDFranLoc_list_log2ChIP[[x]]$sd/sqrt(summaryDFranLoc_list_log2ChIP[[x]]$n-1)
  summaryDFranLoc_list_log2ChIP[[x]]$CI_lower <- summaryDFranLoc_list_log2ChIP[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list_log2ChIP[[x]]$n-1)*summaryDFranLoc_list_log2ChIP[[x]]$sem
  summaryDFranLoc_list_log2ChIP[[x]]$CI_upper <- summaryDFranLoc_list_log2ChIP[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list_log2ChIP[[x]]$n-1)*summaryDFranLoc_list_log2ChIP[[x]]$sem
}

names(summaryDFranLoc_list_log2ChIP) <- ChIPNamesPlot

summaryDFranLoc_list_other  <- mclapply(seq_along(tidyDFranLoc_list_other), function(x) {
  data.frame(window = as.character(wideDFranLoc_list_other[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list_other[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_other[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list_other[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_other[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranLoc_list_other[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_other[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranLoc_list_other))

for(x in seq_along(summaryDFranLoc_list_other)) {
  summaryDFranLoc_list_other[[x]]$window <- factor(summaryDFranLoc_list_other[[x]]$window,
                                                   levels = as.character(wideDFranLoc_list_other[[x]]$window))
  summaryDFranLoc_list_other[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list_other[[x]])[1])
  summaryDFranLoc_list_other[[x]]$sem <- summaryDFranLoc_list_other[[x]]$sd/sqrt(summaryDFranLoc_list_other[[x]]$n-1)
  summaryDFranLoc_list_other[[x]]$CI_lower <- summaryDFranLoc_list_other[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list_other[[x]]$n-1)*summaryDFranLoc_list_other[[x]]$sem
  summaryDFranLoc_list_other[[x]]$CI_upper <- summaryDFranLoc_list_other[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list_other[[x]]$n-1)*summaryDFranLoc_list_other[[x]]$sem
}

names(summaryDFranLoc_list_other) <- otherNamesPlot

summaryDFranLoc_list_control  <- mclapply(seq_along(tidyDFranLoc_list_control), function(x) {
  data.frame(window = as.character(wideDFranLoc_list_control[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list_control[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_control[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list_control[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_control[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranLoc_list_control[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_control[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranLoc_list_control))

for(x in seq_along(summaryDFranLoc_list_control)) {
  summaryDFranLoc_list_control[[x]]$window <- factor(summaryDFranLoc_list_control[[x]]$window,
                                                     levels = as.character(wideDFranLoc_list_control[[x]]$window))
  summaryDFranLoc_list_control[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list_control[[x]])[1])
  summaryDFranLoc_list_control[[x]]$sem <- summaryDFranLoc_list_control[[x]]$sd/sqrt(summaryDFranLoc_list_control[[x]]$n-1)
  summaryDFranLoc_list_control[[x]]$CI_lower <- summaryDFranLoc_list_control[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list_control[[x]]$n-1)*summaryDFranLoc_list_control[[x]]$sem
  summaryDFranLoc_list_control[[x]]$CI_upper <- summaryDFranLoc_list_control[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list_control[[x]]$n-1)*summaryDFranLoc_list_control[[x]]$sem
}

names(summaryDFranLoc_list_control) <- controlNamesPlot

# Convert list summaryDFranLoc_list_log2ChIP into a single data.frame for plotting
summaryDFranLoc_log2ChIP <- bind_rows(summaryDFranLoc_list_log2ChIP, .id = "libName")
summaryDFranLoc_log2ChIP$libName <- factor(summaryDFranLoc_log2ChIP$libName,
                                           levels = names(summaryDFranLoc_list_log2ChIP))

summaryDFranLoc_other <- bind_rows(summaryDFranLoc_list_other, .id = "libName")
summaryDFranLoc_other$libName <- factor(summaryDFranLoc_other$libName,
                                        levels = names(summaryDFranLoc_list_other))

summaryDFranLoc_control <- bind_rows(summaryDFranLoc_list_control, .id = "libName")
summaryDFranLoc_control$libName <- factor(summaryDFranLoc_control$libName,
                                          levels = names(summaryDFranLoc_list_control))


# Define feature start and end labels for plotting
if(grepl("genes", featureName)) {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Define y-axis limits
ymin_list_log2ChIP <- lapply(seq_along(ChIPNamesPlot), function(x) {
  min(c(summaryDFfeature_log2ChIP[summaryDFfeature_log2ChIP$libName ==
                                    ChIPNamesPlot[x],]$CI_lower,
        summaryDFranLoc_log2ChIP[summaryDFranLoc_log2ChIP$libName ==
                                   ChIPNamesPlot[x],]$CI_lower))
})
ymax_list_log2ChIP <- lapply(seq_along(ChIPNamesPlot), function(x) {
  max(c(summaryDFfeature_log2ChIP[summaryDFfeature_log2ChIP$libName ==
                                    ChIPNamesPlot[x],]$CI_upper,
        summaryDFranLoc_log2ChIP[summaryDFranLoc_log2ChIP$libName ==
                                   ChIPNamesPlot[x],]$CI_upper))
})

ymin_list_other <- lapply(seq_along(otherNamesPlot), function(x) {
  min(c(summaryDFfeature_other[summaryDFfeature_other$libName ==
                                 otherNamesPlot[x],]$CI_lower,
        summaryDFranLoc_other[summaryDFranLoc_other$libName ==
                                otherNamesPlot[x],]$CI_lower))
})
ymax_list_other <- lapply(seq_along(otherNamesPlot), function(x) {
  max(c(summaryDFfeature_other[summaryDFfeature_other$libName ==
                                 otherNamesPlot[x],]$CI_upper,
        summaryDFranLoc_other[summaryDFranLoc_other$libName ==
                                otherNamesPlot[x],]$CI_upper))
})

ymin_list_control <- lapply(seq_along(controlNamesPlot), function(x) {
  min(c(summaryDFfeature_control[summaryDFfeature_control$libName ==
                                   controlNamesPlot[x],]$CI_lower,
        summaryDFranLoc_control[summaryDFranLoc_control$libName ==
                                  controlNamesPlot[x],]$CI_lower))
})
ymax_list_control <- lapply(seq_along(controlNamesPlot), function(x) {
  max(c(summaryDFfeature_control[summaryDFfeature_control$libName ==
                                   controlNamesPlot[x],]$CI_upper,
        summaryDFranLoc_control[summaryDFranLoc_control$libName ==
                                  controlNamesPlot[x],]$CI_upper))
})

# Plot average coverage profiles with 95% CI ribbon
## feature
ggObjGA_combined_log2ChIP <- NULL
ggObj1_combined_log2ChIP <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[summaryDFfeature_log2ChIP$libName ==
                                                  ChIPNamesPlot[x],]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = ChIPColours[x]) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = ChIPColours[x]) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_list_log2ChIP[[x]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_list_log2ChIP[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_list_log2ChIP[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.1,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranLoc
ggObj2_combined_log2ChIP <- lapply(seq_along(summaryDFranLoc_list_log2ChIP), function(x) {
  summaryDFranLoc <- summaryDFranLoc_log2ChIP[summaryDFranLoc_log2ChIP$libName ==
                                                  ChIPNamesPlot[x],]
  ggplot(data = summaryDFranLoc,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranLoc,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = ChIPColours[x]) +
  geom_ribbon(data = summaryDFranLoc,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = ChIPColours[x]) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranLoc_list_log2ChIP[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranLoc_list_log2ChIP[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranLoc_list_log2ChIP[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.1,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("Random loci (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranLoc$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
ggObjGA_combined_log2ChIP <- grid.arrange(grobs = c(ggObj1_combined_log2ChIP,
                                                    ggObj2_combined_log2ChIP),
                                          layout_matrix = cbind(1:length(ChIPNames),
                                                                (length(ChIPNames)+1):(length(ChIPNames)*2)))
                                          #nrow = length(ChIPNames), ncol = 2)
ggsave(paste0(plotDir,
              "avgProfiles_around_",
              featureNamePlot, "_in_cluster", clusterNo,
              "_by_log2_ASY1_CS_Rep1_ChIP_control_in_", region,
              "_of_", featureName, ".pdf"),
       plot = ggObjGA_combined_log2ChIP,
       height = 6.5*length(ChIPNames), width = 14, limitsize = FALSE)
