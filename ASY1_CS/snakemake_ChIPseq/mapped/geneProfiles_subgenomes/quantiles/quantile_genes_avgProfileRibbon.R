#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage profiles with 95% CIs around
# gene quantiles; e.g.,
# clusters_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters/cluster1_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters_of_genes_in_Agenome_genomewide.txt

# Usage:
# /applications/R/R-3.5.0/bin/Rscript quantile_genes_avgProfileRibbon.R ASY1_CS_Rep1_ChIP ASY1_CS both 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' 3500 2000 2kb '2 kb' 20 20bp promoters 4 100kb 1

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#align <- "both"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#bodyLength <- 3500
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "promoters"
#quantiles <- 4
#winName <- "100kb"
#minMarkerDist <- 1

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

outDir <- paste0("quantiles_by_log2_", libName,
                 "_control_in_", region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]

# Load table of features grouped into quantiles
# by decreasing log2(libName/control)
featuresDF <- read.table(paste0(outDir, "features_", quantiles, "quantiles",
                                "_by_log2_", libName, "_control_in_",
                                region, "_of_",
                                substring(featureName[1][1], first = 1, last = 5), "_in_",
                                paste0(substring(featureName, first = 10, last = 16),
                                       collapse = "_"), "_",
                                substring(featureName[1][1], first = 18), ".txt"),
                         header = T, sep = "\t")

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
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
stopifnot(identical(as.character(featuresDF$featureID),
                    as.character(features$V9)))

# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$quantile == paste0("Quantile ", k))
})

## Random feature quantiles
# Define function to randomly select n rows from
# a data.frame
selectRandomFeatures <- function(features, n) {
  return(features[sample(x = dim(features)[1],
                         size = n,
                         replace = FALSE),])
}

# Define seed so that random selections are reproducible
set.seed(93750174)

# Divide features into random sets of equal number,
# with the same number of genes per chromosome as
# above-defined libName-defined feature quantiles
randomPCIndices <- lapply(1:quantiles, function(k) {
  randomPCIndicesk <- NULL
  for(i in 1:length(chrs)) {
    randomPCfeatureskChr <- selectRandomFeatures(features = featuresDF[featuresDF$seqnames == chrs[i],],
                                                 n = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k) &
                                                                    featuresDF$seqnames == chrs[i],])[1])
    randomPCIndicesk <- c(randomPCIndicesk, as.integer(rownames(randomPCfeatureskChr)))
  }
  randomPCIndicesk
})
# Confirm per-chromosome feature numbers are the same for quantiles and random groupings
sapply(seq_along(chrs), function(x) {
  if(!identical(dim(featuresDF[randomPCIndices[[1]],][featuresDF$seqnames == chrs[1],]),
                dim(featuresDF[quantileIndices[[1]],][featuresDF$seqnames == chrs[1],]))) {
    stop("Quantile features and random features do not consist of the same number of features per chromosome")
  }
})


## Load feature matrices (featureMats) for each control library
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
                                featureName[y], "_matrix_bin", binName,
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
                "DNaseI_Rep1_SRR8447247",
                "WT_RNAseq_Rep1_ERR2402974",
                "WT_RNAseq_Rep2_ERR2402973",
                "WT_RNAseq_Rep3_ERR2402972"
               )
otherNamesDir <- c(
                   "MNase",
                   "DNaseI",
                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci"
                  )
otherNamesPlot <- c(
                    "MNase",
                    "DNaseI",
                    "RNA-seq Rep1",
                    "RNA-seq Rep2",
                    "RNA-seq Rep3"
                   )
otherColours <- c(
                  "darkcyan",
                  "purple",
                  "red4",
                  "red4",
                  "red4"
                 )
sRNANames <- c(
               "CS+_2_LIB18613_LDI16228"
              )
sRNANamesDir <- c(
                  "sRNAseq_meiocyte_Martin_Moore"
                 )
sRNANamesPlot <- c(
                   "20-nt sRNAs",
                   "21-nt sRNAs",
                   "22-nt sRNAs",
                   "23-nt sRNAs",
                   "24-nt sRNAs",
                   "34-nt sRNAs"
                  )
sRNAsizes <- c(
               "20nt",
               "21nt",
               "22nt",
               "23nt",
               "24nt",
               "33nt",
               "34nt"
              )
sRNAColours <- c(
                 "red",
                 "blue",
                 "green2",
                 "darkorange2",
                 "purple3",
                 "darkgreen",
                 "deeppink"
                )
DNAmethNames <- c(
                  "BSseq_Rep8a_SRR6792678"
                 )
DNAmethNamesDir <- c(
                     "BSseq"
                    )
DNAmethContexts <- c(
                     "CpG",
                     "CHG",
                     "CHH"
                    )
DNAmethNamesPlot <- c(
                      "mCG",
                      "mCHG",
                      "mCHH"
                     )
DNAmethColours <- c(
                    "navy",
                    "blue",
                    "deepskyblue1"
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
  } else if(grepl("RNAseq", otherNames[x])) {
    paste0("/home/ajt200/analysis/wheat/",
           otherNamesDir[x], "/snakemake_RNAseq_HISAT2/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("otherNames[", x, "] is not compatible with the specified coverage matrix paths"))
  }
})
sRNADirs <- sapply(seq_along(sRNANames), function(x) {
  if(sRNANames[x] %in% c("CS+_2_LIB18613_LDI16228")) {
    paste0("/home/ajt200/analysis/wheat/",
           sRNANamesDir[x], "/snakemake_sRNAseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("sRNANames[", x, "] is not compatible with the specified coverage matrix paths"))
  }
})
DNAmethDirs <- sapply(seq_along(DNAmethNames), function(x) {
  if(DNAmethNames[x] %in% c("BSseq_Rep8a_SRR6792678")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           DNAmethNamesDir[x],
           "/snakemake_BSseq/coverage/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("DNAmethNames[", x, "] is not compatible with the specified coverage matrix paths"))
  }
})

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

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
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

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
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


# Add column names
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
log2ChIP_mats_quantiles <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  list(
       # feature
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][quantileIndices[[k]],]
       }),
       # ranFeat
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][randomPCIndices[[k]],]
       }),
       # ranLoc
       lapply(1:quantiles, function(k) {
         log2ChIP_ranLocMats[[x]][quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(log2ChIP_featureMats))

## feature
# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats_quantiles), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(log2ChIP_mats_quantiles[[x]][[y]][[k]]),
                 t(log2ChIP_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(log2ChIP_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
    })
  })
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]])[1])
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
    })
  })
}

names(summaryDFfeature_list_log2ChIP) <- ChIPNamesPlot

# Convert list summaryDFfeature_list_log2ChIP into a single data.frame for plotting
summaryDFfeature_log2ChIP <- bind_rows(summaryDFfeature_list_log2ChIP, .id = "libName")
summaryDFfeature_log2ChIP$libName <- factor(summaryDFfeature_log2ChIP$libName,
                                            levels = names(summaryDFfeature_list_log2ChIP))

summaryDFfeature_other <- bind_rows(summaryDFfeature_list_other, .id = "libName")
summaryDFfeature_other$libName <- factor(summaryDFfeature_other$libName,
                                         levels = names(summaryDFfeature_list_other))

summaryDFfeature_sRNA <- bind_rows(summaryDFfeature_list_sRNA, .id = "libName")
summaryDFfeature_sRNA$libName <- factor(summaryDFfeature_sRNA$libName,
                                        levels = names(summaryDFfeature_list_sRNA))

summaryDFfeature_control <- bind_rows(summaryDFfeature_list_control, .id = "libName")
summaryDFfeature_control$libName <- factor(summaryDFfeature_control$libName,
                                           levels = names(summaryDFfeature_list_control))

summaryDFfeature_DNAmeth <- bind_rows(summaryDFfeature_list_DNAmeth, .id = "libName")
summaryDFfeature_DNAmeth$libName <- factor(summaryDFfeature_DNAmeth$libName,
                                           levels = names(summaryDFfeature_list_DNAmeth))

## ranFeat
# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFranFeat_list_log2ChIP <- mclapply(seq_along(log2ChIP_ranFeatMats), function(x) {
  data.frame(window = colnames(log2ChIP_ranFeatMats[[x]]),
             t(log2ChIP_ranFeatMats[[x]]))
}, mc.cores = length(log2ChIP_ranFeatMats))

wideDFranFeat_list_other <- mclapply(seq_along(other_ranFeatMats), function(x) {
  data.frame(window = colnames(other_ranFeatMats[[x]]),
             t(other_ranFeatMats[[x]]))
}, mc.cores = length(other_ranFeatMats))

wideDFranFeat_list_sRNA <- mclapply(seq_along(sRNA_ranFeatMats), function(x) {
  data.frame(window = colnames(sRNA_ranFeatMats[[x]]),
             t(sRNA_ranFeatMats[[x]]))
}, mc.cores = length(sRNA_ranFeatMats))

wideDFranFeat_list_control <- mclapply(seq_along(control_ranFeatMats), function(x) {
  data.frame(window = colnames(control_ranFeatMats[[x]]),
             t(control_ranFeatMats[[x]]))
}, mc.cores = length(control_ranFeatMats))

wideDFranFeat_list_DNAmeth <- mclapply(seq_along(DNAmeth_ranFeatMats), function(x) {
  data.frame(window = colnames(DNAmeth_ranFeatMats[[x]]),
             t(DNAmeth_ranFeatMats[[x]]))
}, mc.cores = length(DNAmeth_ranFeatMats))

# Convert into tidy data.frame (long format)
tidyDFranFeat_list_log2ChIP  <- mclapply(seq_along(wideDFranFeat_list_log2ChIP), function(x) {
  gather(data  = wideDFranFeat_list_log2ChIP[[x]],
         key   = ranFeat,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranFeat_list_log2ChIP))

tidyDFranFeat_list_other  <- mclapply(seq_along(wideDFranFeat_list_other), function(x) {
  gather(data  = wideDFranFeat_list_other[[x]],
         key   = ranFeat,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranFeat_list_other))

tidyDFranFeat_list_sRNA  <- mclapply(seq_along(wideDFranFeat_list_sRNA), function(x) {
  gather(data  = wideDFranFeat_list_sRNA[[x]],
         key   = ranFeat,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranFeat_list_sRNA))

tidyDFranFeat_list_control  <- mclapply(seq_along(wideDFranFeat_list_control), function(x) {
  gather(data  = wideDFranFeat_list_control[[x]],
         key   = ranFeat,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranFeat_list_control))

tidyDFranFeat_list_DNAmeth  <- mclapply(seq_along(wideDFranFeat_list_DNAmeth), function(x) {
  gather(data  = wideDFranFeat_list_DNAmeth[[x]],
         key   = ranFeat,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranFeat_list_DNAmeth))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFranFeat_list_log2ChIP)) {
  tidyDFranFeat_list_log2ChIP[[x]]$window <- factor(tidyDFranFeat_list_log2ChIP[[x]]$window,
                                                    levels = as.character(wideDFranFeat_list_log2ChIP[[x]]$window))
}

for(x in seq_along(tidyDFranFeat_list_other)) {
  tidyDFranFeat_list_other[[x]]$window <- factor(tidyDFranFeat_list_other[[x]]$window,
                                                 levels = as.character(wideDFranFeat_list_other[[x]]$window))
}

for(x in seq_along(tidyDFranFeat_list_sRNA)) {
  tidyDFranFeat_list_sRNA[[x]]$window <- factor(tidyDFranFeat_list_sRNA[[x]]$window,
                                                levels = as.character(wideDFranFeat_list_sRNA[[x]]$window))
}

for(x in seq_along(tidyDFranFeat_list_control)) {
  tidyDFranFeat_list_control[[x]]$window <- factor(tidyDFranFeat_list_control[[x]]$window,
                                                   levels = as.character(wideDFranFeat_list_control[[x]]$window))
}

for(x in seq_along(tidyDFranFeat_list_DNAmeth)) {
  tidyDFranFeat_list_DNAmeth[[x]]$window <- factor(tidyDFranFeat_list_DNAmeth[[x]]$window,
                                                   levels = as.character(wideDFranFeat_list_DNAmeth[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (ranFeats) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFranFeat_list_log2ChIP  <- mclapply(seq_along(tidyDFranFeat_list_log2ChIP), function(x) {
  data.frame(window = as.character(wideDFranFeat_list_log2ChIP[[x]]$window),
             n      = tapply(X     = tidyDFranFeat_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_log2ChIP[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranFeat_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_log2ChIP[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranFeat_list_log2ChIP[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_log2ChIP[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranFeat_list_log2ChIP))

for(x in seq_along(summaryDFranFeat_list_log2ChIP)) {
  summaryDFranFeat_list_log2ChIP[[x]]$window <- factor(summaryDFranFeat_list_log2ChIP[[x]]$window,
                                                       levels = as.character(wideDFranFeat_list_log2ChIP[[x]]$window))
  summaryDFranFeat_list_log2ChIP[[x]]$winNo <- factor(1:dim(summaryDFranFeat_list_log2ChIP[[x]])[1])
  summaryDFranFeat_list_log2ChIP[[x]]$sem <- summaryDFranFeat_list_log2ChIP[[x]]$sd/sqrt(summaryDFranFeat_list_log2ChIP[[x]]$n-1)
  summaryDFranFeat_list_log2ChIP[[x]]$CI_lower <- summaryDFranFeat_list_log2ChIP[[x]]$mean -
    qt(0.975, df = summaryDFranFeat_list_log2ChIP[[x]]$n-1)*summaryDFranFeat_list_log2ChIP[[x]]$sem
  summaryDFranFeat_list_log2ChIP[[x]]$CI_upper <- summaryDFranFeat_list_log2ChIP[[x]]$mean +
    qt(0.975, df = summaryDFranFeat_list_log2ChIP[[x]]$n-1)*summaryDFranFeat_list_log2ChIP[[x]]$sem
}

names(summaryDFranFeat_list_log2ChIP) <- ChIPNamesPlot

summaryDFranFeat_list_other  <- mclapply(seq_along(tidyDFranFeat_list_other), function(x) {
  data.frame(window = as.character(wideDFranFeat_list_other[[x]]$window),
             n      = tapply(X     = tidyDFranFeat_list_other[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_other[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranFeat_list_other[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_other[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranFeat_list_other[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_other[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranFeat_list_other))

for(x in seq_along(summaryDFranFeat_list_other)) {
  summaryDFranFeat_list_other[[x]]$window <- factor(summaryDFranFeat_list_other[[x]]$window,
                                                    levels = as.character(wideDFranFeat_list_other[[x]]$window))
  summaryDFranFeat_list_other[[x]]$winNo <- factor(1:dim(summaryDFranFeat_list_other[[x]])[1])
  summaryDFranFeat_list_other[[x]]$sem <- summaryDFranFeat_list_other[[x]]$sd/sqrt(summaryDFranFeat_list_other[[x]]$n-1)
  summaryDFranFeat_list_other[[x]]$CI_lower <- summaryDFranFeat_list_other[[x]]$mean -
    qt(0.975, df = summaryDFranFeat_list_other[[x]]$n-1)*summaryDFranFeat_list_other[[x]]$sem
  summaryDFranFeat_list_other[[x]]$CI_upper <- summaryDFranFeat_list_other[[x]]$mean +
    qt(0.975, df = summaryDFranFeat_list_other[[x]]$n-1)*summaryDFranFeat_list_other[[x]]$sem
}

names(summaryDFranFeat_list_other) <- otherNamesPlot

summaryDFranFeat_list_sRNA  <- mclapply(seq_along(tidyDFranFeat_list_sRNA), function(x) {
  data.frame(window = as.character(wideDFranFeat_list_sRNA[[x]]$window),
             n      = tapply(X     = tidyDFranFeat_list_sRNA[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_sRNA[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranFeat_list_sRNA[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_sRNA[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranFeat_list_sRNA[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_sRNA[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranFeat_list_sRNA))

for(x in seq_along(summaryDFranFeat_list_sRNA)) {
  summaryDFranFeat_list_sRNA[[x]]$window <- factor(summaryDFranFeat_list_sRNA[[x]]$window,
                                                    levels = as.character(wideDFranFeat_list_sRNA[[x]]$window))
  summaryDFranFeat_list_sRNA[[x]]$winNo <- factor(1:dim(summaryDFranFeat_list_sRNA[[x]])[1])
  summaryDFranFeat_list_sRNA[[x]]$sem <- summaryDFranFeat_list_sRNA[[x]]$sd/sqrt(summaryDFranFeat_list_sRNA[[x]]$n-1)
  summaryDFranFeat_list_sRNA[[x]]$CI_lower <- summaryDFranFeat_list_sRNA[[x]]$mean -
    qt(0.975, df = summaryDFranFeat_list_sRNA[[x]]$n-1)*summaryDFranFeat_list_sRNA[[x]]$sem
  summaryDFranFeat_list_sRNA[[x]]$CI_upper <- summaryDFranFeat_list_sRNA[[x]]$mean +
    qt(0.975, df = summaryDFranFeat_list_sRNA[[x]]$n-1)*summaryDFranFeat_list_sRNA[[x]]$sem
}

names(summaryDFranFeat_list_sRNA) <- sRNANamesPlot

summaryDFranFeat_list_control  <- mclapply(seq_along(tidyDFranFeat_list_control), function(x) {
  data.frame(window = as.character(wideDFranFeat_list_control[[x]]$window),
             n      = tapply(X     = tidyDFranFeat_list_control[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_control[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranFeat_list_control[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_control[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranFeat_list_control[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_control[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranFeat_list_control))

for(x in seq_along(summaryDFranFeat_list_control)) {
  summaryDFranFeat_list_control[[x]]$window <- factor(summaryDFranFeat_list_control[[x]]$window,
                                                      levels = as.character(wideDFranFeat_list_control[[x]]$window))
  summaryDFranFeat_list_control[[x]]$winNo <- factor(1:dim(summaryDFranFeat_list_control[[x]])[1])
  summaryDFranFeat_list_control[[x]]$sem <- summaryDFranFeat_list_control[[x]]$sd/sqrt(summaryDFranFeat_list_control[[x]]$n-1)
  summaryDFranFeat_list_control[[x]]$CI_lower <- summaryDFranFeat_list_control[[x]]$mean -
    qt(0.975, df = summaryDFranFeat_list_control[[x]]$n-1)*summaryDFranFeat_list_control[[x]]$sem
  summaryDFranFeat_list_control[[x]]$CI_upper <- summaryDFranFeat_list_control[[x]]$mean +
    qt(0.975, df = summaryDFranFeat_list_control[[x]]$n-1)*summaryDFranFeat_list_control[[x]]$sem
}

names(summaryDFranFeat_list_control) <- controlNamesPlot

summaryDFranFeat_list_DNAmeth  <- mclapply(seq_along(tidyDFranFeat_list_DNAmeth), function(x) {
  data.frame(window = as.character(wideDFranFeat_list_DNAmeth[[x]]$window),
             n      = tapply(X     = tidyDFranFeat_list_DNAmeth[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_DNAmeth[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranFeat_list_DNAmeth[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_DNAmeth[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranFeat_list_DNAmeth[[x]]$coverage,
                             INDEX = tidyDFranFeat_list_DNAmeth[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranFeat_list_DNAmeth))

for(x in seq_along(summaryDFranFeat_list_DNAmeth)) {
  summaryDFranFeat_list_DNAmeth[[x]]$window <- factor(summaryDFranFeat_list_DNAmeth[[x]]$window,
                                                      levels = as.character(wideDFranFeat_list_DNAmeth[[x]]$window))
  summaryDFranFeat_list_DNAmeth[[x]]$winNo <- factor(1:dim(summaryDFranFeat_list_DNAmeth[[x]])[1])
  summaryDFranFeat_list_DNAmeth[[x]]$sem <- summaryDFranFeat_list_DNAmeth[[x]]$sd/sqrt(summaryDFranFeat_list_DNAmeth[[x]]$n-1)
  summaryDFranFeat_list_DNAmeth[[x]]$CI_lower <- summaryDFranFeat_list_DNAmeth[[x]]$mean -
    qt(0.975, df = summaryDFranFeat_list_DNAmeth[[x]]$n-1)*summaryDFranFeat_list_DNAmeth[[x]]$sem
  summaryDFranFeat_list_DNAmeth[[x]]$CI_upper <- summaryDFranFeat_list_DNAmeth[[x]]$mean +
    qt(0.975, df = summaryDFranFeat_list_DNAmeth[[x]]$n-1)*summaryDFranFeat_list_DNAmeth[[x]]$sem
}

names(summaryDFranFeat_list_DNAmeth) <- DNAmethNamesPlot

# Convert list summaryDFranFeat_list_log2ChIP into a single data.frame for plotting
summaryDFranFeat_log2ChIP <- bind_rows(summaryDFranFeat_list_log2ChIP, .id = "libName")
summaryDFranFeat_log2ChIP$libName <- factor(summaryDFranFeat_log2ChIP$libName,
                                            levels = names(summaryDFranFeat_list_log2ChIP))

summaryDFranFeat_other <- bind_rows(summaryDFranFeat_list_other, .id = "libName")
summaryDFranFeat_other$libName <- factor(summaryDFranFeat_other$libName,
                                         levels = names(summaryDFranFeat_list_other))

summaryDFranFeat_sRNA <- bind_rows(summaryDFranFeat_list_sRNA, .id = "libName")
summaryDFranFeat_sRNA$libName <- factor(summaryDFranFeat_sRNA$libName,
                                        levels = names(summaryDFranFeat_list_sRNA))

summaryDFranFeat_control <- bind_rows(summaryDFranFeat_list_control, .id = "libName")
summaryDFranFeat_control$libName <- factor(summaryDFranFeat_control$libName,
                                           levels = names(summaryDFranFeat_list_control))

summaryDFranFeat_DNAmeth <- bind_rows(summaryDFranFeat_list_DNAmeth, .id = "libName")
summaryDFranFeat_DNAmeth$libName <- factor(summaryDFranFeat_DNAmeth$libName,
                                           levels = names(summaryDFranFeat_list_DNAmeth))

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

wideDFranLoc_list_sRNA <- mclapply(seq_along(sRNA_ranLocMats), function(x) {
  data.frame(window = colnames(sRNA_ranLocMats[[x]]),
             t(sRNA_ranLocMats[[x]]))
}, mc.cores = length(sRNA_ranLocMats))

wideDFranLoc_list_control <- mclapply(seq_along(control_ranLocMats), function(x) {
  data.frame(window = colnames(control_ranLocMats[[x]]),
             t(control_ranLocMats[[x]]))
}, mc.cores = length(control_ranLocMats))

wideDFranLoc_list_DNAmeth <- mclapply(seq_along(DNAmeth_ranLocMats), function(x) {
  data.frame(window = colnames(DNAmeth_ranLocMats[[x]]),
             t(DNAmeth_ranLocMats[[x]]))
}, mc.cores = length(DNAmeth_ranLocMats))

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

tidyDFranLoc_list_sRNA  <- mclapply(seq_along(wideDFranLoc_list_sRNA), function(x) {
  gather(data  = wideDFranLoc_list_sRNA[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list_sRNA))

tidyDFranLoc_list_control  <- mclapply(seq_along(wideDFranLoc_list_control), function(x) {
  gather(data  = wideDFranLoc_list_control[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list_control))

tidyDFranLoc_list_DNAmeth  <- mclapply(seq_along(wideDFranLoc_list_DNAmeth), function(x) {
  gather(data  = wideDFranLoc_list_DNAmeth[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list_DNAmeth))

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

for(x in seq_along(tidyDFranLoc_list_sRNA)) {
  tidyDFranLoc_list_sRNA[[x]]$window <- factor(tidyDFranLoc_list_sRNA[[x]]$window,
                                               levels = as.character(wideDFranLoc_list_sRNA[[x]]$window))
}

for(x in seq_along(tidyDFranLoc_list_control)) {
  tidyDFranLoc_list_control[[x]]$window <- factor(tidyDFranLoc_list_control[[x]]$window,
                                                  levels = as.character(wideDFranLoc_list_control[[x]]$window))
}

for(x in seq_along(tidyDFranLoc_list_DNAmeth)) {
  tidyDFranLoc_list_DNAmeth[[x]]$window <- factor(tidyDFranLoc_list_DNAmeth[[x]]$window,
                                                  levels = as.character(wideDFranLoc_list_DNAmeth[[x]]$window))
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

summaryDFranLoc_list_sRNA  <- mclapply(seq_along(tidyDFranLoc_list_sRNA), function(x) {
  data.frame(window = as.character(wideDFranLoc_list_sRNA[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list_sRNA[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_sRNA[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list_sRNA[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_sRNA[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranLoc_list_sRNA[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_sRNA[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranLoc_list_sRNA))

for(x in seq_along(summaryDFranLoc_list_sRNA)) {
  summaryDFranLoc_list_sRNA[[x]]$window <- factor(summaryDFranLoc_list_sRNA[[x]]$window,
                                                   levels = as.character(wideDFranLoc_list_sRNA[[x]]$window))
  summaryDFranLoc_list_sRNA[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list_sRNA[[x]])[1])
  summaryDFranLoc_list_sRNA[[x]]$sem <- summaryDFranLoc_list_sRNA[[x]]$sd/sqrt(summaryDFranLoc_list_sRNA[[x]]$n-1)
  summaryDFranLoc_list_sRNA[[x]]$CI_lower <- summaryDFranLoc_list_sRNA[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list_sRNA[[x]]$n-1)*summaryDFranLoc_list_sRNA[[x]]$sem
  summaryDFranLoc_list_sRNA[[x]]$CI_upper <- summaryDFranLoc_list_sRNA[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list_sRNA[[x]]$n-1)*summaryDFranLoc_list_sRNA[[x]]$sem
}

names(summaryDFranLoc_list_sRNA) <- sRNANamesPlot

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

summaryDFranLoc_list_DNAmeth  <- mclapply(seq_along(tidyDFranLoc_list_DNAmeth), function(x) {
  data.frame(window = as.character(wideDFranLoc_list_DNAmeth[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list_DNAmeth[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_DNAmeth[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list_DNAmeth[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_DNAmeth[[x]]$window,
                             FUN   = mean,
                             na.rm = TRUE),
             sd     = tapply(X     = tidyDFranLoc_list_DNAmeth[[x]]$coverage,
                             INDEX = tidyDFranLoc_list_DNAmeth[[x]]$window,
                             FUN   = sd,
                             na.rm = TRUE))
}, mc.cores = length(tidyDFranLoc_list_DNAmeth))

for(x in seq_along(summaryDFranLoc_list_DNAmeth)) {
  summaryDFranLoc_list_DNAmeth[[x]]$window <- factor(summaryDFranLoc_list_DNAmeth[[x]]$window,
                                                     levels = as.character(wideDFranLoc_list_DNAmeth[[x]]$window))
  summaryDFranLoc_list_DNAmeth[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list_DNAmeth[[x]])[1])
  summaryDFranLoc_list_DNAmeth[[x]]$sem <- summaryDFranLoc_list_DNAmeth[[x]]$sd/sqrt(summaryDFranLoc_list_DNAmeth[[x]]$n-1)
  summaryDFranLoc_list_DNAmeth[[x]]$CI_lower <- summaryDFranLoc_list_DNAmeth[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list_DNAmeth[[x]]$n-1)*summaryDFranLoc_list_DNAmeth[[x]]$sem
  summaryDFranLoc_list_DNAmeth[[x]]$CI_upper <- summaryDFranLoc_list_DNAmeth[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list_DNAmeth[[x]]$n-1)*summaryDFranLoc_list_DNAmeth[[x]]$sem
}

names(summaryDFranLoc_list_DNAmeth) <- DNAmethNamesPlot

# Convert list summaryDFranLoc_list_log2ChIP into a single data.frame for plotting
summaryDFranLoc_log2ChIP <- bind_rows(summaryDFranLoc_list_log2ChIP, .id = "libName")
summaryDFranLoc_log2ChIP$libName <- factor(summaryDFranLoc_log2ChIP$libName,
                                           levels = names(summaryDFranLoc_list_log2ChIP))

summaryDFranLoc_other <- bind_rows(summaryDFranLoc_list_other, .id = "libName")
summaryDFranLoc_other$libName <- factor(summaryDFranLoc_other$libName,
                                        levels = names(summaryDFranLoc_list_other))

summaryDFranLoc_sRNA <- bind_rows(summaryDFranLoc_list_sRNA, .id = "libName")
summaryDFranLoc_sRNA$libName <- factor(summaryDFranLoc_sRNA$libName,
                                       levels = names(summaryDFranLoc_list_sRNA))

summaryDFranLoc_control <- bind_rows(summaryDFranLoc_list_control, .id = "libName")
summaryDFranLoc_control$libName <- factor(summaryDFranLoc_control$libName,
                                          levels = names(summaryDFranLoc_list_control))

summaryDFranLoc_DNAmeth <- bind_rows(summaryDFranLoc_list_DNAmeth, .id = "libName")
summaryDFranLoc_DNAmeth$libName <- factor(summaryDFranLoc_DNAmeth$libName,
                                          levels = names(summaryDFranLoc_list_DNAmeth))


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
        summaryDFranFeat_log2ChIP[summaryDFranFeat_log2ChIP$libName ==
                                    ChIPNamesPlot[x],]$CI_lower,
        summaryDFranLoc_log2ChIP[summaryDFranLoc_log2ChIP$libName ==
                                   ChIPNamesPlot[x],]$CI_lower))
})
ymax_list_log2ChIP <- lapply(seq_along(ChIPNamesPlot), function(x) {
  max(c(summaryDFfeature_log2ChIP[summaryDFfeature_log2ChIP$libName ==
                                    ChIPNamesPlot[x],]$CI_upper,
        summaryDFranFeat_log2ChIP[summaryDFranFeat_log2ChIP$libName ==
                                    ChIPNamesPlot[x],]$CI_upper,
        summaryDFranLoc_log2ChIP[summaryDFranLoc_log2ChIP$libName ==
                                   ChIPNamesPlot[x],]$CI_upper))
})

ymin_list_other <- lapply(seq_along(otherNamesPlot), function(x) {
  min(c(summaryDFfeature_other[summaryDFfeature_other$libName ==
                                 otherNamesPlot[x],]$CI_lower,
        summaryDFranFeat_other[summaryDFranFeat_other$libName ==
                                 otherNamesPlot[x],]$CI_lower,
        summaryDFranLoc_other[summaryDFranLoc_other$libName ==
                                otherNamesPlot[x],]$CI_lower))
})
ymax_list_other <- lapply(seq_along(otherNamesPlot), function(x) {
  max(c(summaryDFfeature_other[summaryDFfeature_other$libName ==
                                 otherNamesPlot[x],]$CI_upper,
        summaryDFranFeat_other[summaryDFranFeat_other$libName ==
                                 otherNamesPlot[x],]$CI_upper,
        summaryDFranLoc_other[summaryDFranLoc_other$libName ==
                                otherNamesPlot[x],]$CI_upper))
})

ymin_list_sRNA <- lapply(seq_along(sRNANamesPlot), function(x) {
  min(c(summaryDFfeature_sRNA[summaryDFfeature_sRNA$libName ==
                                sRNANamesPlot[x],]$CI_lower,
        summaryDFranFeat_sRNA[summaryDFranFeat_sRNA$libName ==
                                sRNANamesPlot[x],]$CI_lower,
        summaryDFranLoc_sRNA[summaryDFranLoc_sRNA$libName ==
                               sRNANamesPlot[x],]$CI_lower))
})
ymax_list_sRNA <- lapply(seq_along(sRNANamesPlot), function(x) {
  max(c(summaryDFfeature_sRNA[summaryDFfeature_sRNA$libName ==
                                sRNANamesPlot[x],]$CI_upper,
        summaryDFranFeat_sRNA[summaryDFranFeat_sRNA$libName ==
                                sRNANamesPlot[x],]$CI_upper,
        summaryDFranLoc_sRNA[summaryDFranLoc_sRNA$libName ==
                               sRNANamesPlot[x],]$CI_upper))
})

ymin_list_control <- lapply(seq_along(controlNamesPlot), function(x) {
  min(c(summaryDFfeature_control[summaryDFfeature_control$libName ==
                                   controlNamesPlot[x],]$CI_lower,
        summaryDFranFeat_control[summaryDFranFeat_control$libName ==
                                   controlNamesPlot[x],]$CI_lower,
        summaryDFranLoc_control[summaryDFranLoc_control$libName ==
                                  controlNamesPlot[x],]$CI_lower))
})
ymax_list_control <- lapply(seq_along(controlNamesPlot), function(x) {
  max(c(summaryDFfeature_control[summaryDFfeature_control$libName ==
                                   controlNamesPlot[x],]$CI_upper,
        summaryDFranFeat_control[summaryDFranFeat_control$libName ==
                                   controlNamesPlot[x],]$CI_upper,
        summaryDFranLoc_control[summaryDFranLoc_control$libName ==
                                  controlNamesPlot[x],]$CI_upper))
})

ymin_list_DNAmeth <- lapply(seq_along(DNAmethNamesPlot), function(x) {
  min(c(summaryDFfeature_DNAmeth[summaryDFfeature_DNAmeth$libName ==
                                   DNAmethNamesPlot[x],]$CI_lower,
        summaryDFranFeat_DNAmeth[summaryDFranFeat_DNAmeth$libName ==
                                   DNAmethNamesPlot[x],]$CI_lower,
        summaryDFranLoc_DNAmeth[summaryDFranLoc_DNAmeth$libName ==
                                  DNAmethNamesPlot[x],]$CI_lower))
})
ymax_list_DNAmeth <- lapply(seq_along(DNAmethNamesPlot), function(x) {
  max(c(summaryDFfeature_DNAmeth[summaryDFfeature_DNAmeth$libName ==
                                   DNAmethNamesPlot[x],]$CI_upper,
        summaryDFranFeat_DNAmeth[summaryDFranFeat_DNAmeth$libName ==
                                   DNAmethNamesPlot[x],]$CI_upper,
        summaryDFranLoc_DNAmeth[summaryDFranLoc_DNAmeth$libName ==
                                  DNAmethNamesPlot[x],]$CI_upper))
})

# Plot average coverage profiles with 95% CI ribbon
## feature
ggObjGA_combined <- NULL
ggObj1_combined_log2ChIP <- lapply(seq_along(ChIPNames), function(x) {
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
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranFeat
ggObj2_combined_log2ChIP <- lapply(seq_along(ChIPNames), function(x) {
  summaryDFranFeat <- summaryDFranFeat_log2ChIP[summaryDFranFeat_log2ChIP$libName ==
                                                  ChIPNamesPlot[x],]
  ggplot(data = summaryDFranFeat,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranFeat,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = ChIPColours[x]) +
  geom_ribbon(data = summaryDFranFeat,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = ChIPColours[x]) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranFeat_list_log2ChIP[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranFeat_list_log2ChIP[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranFeat_list_log2ChIP[[x]])[1])-(downstream/binSize)),
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
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("ASY1 Cluster" * .(clusterLast) ~ "genes (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranFeat$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranLoc
ggObj3_combined_log2ChIP <- lapply(seq_along(ChIPNames), function(x) {
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
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("Random loci (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranLoc$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## feature
ggObj1_combined_other <- lapply(seq_along(otherNames), function(x) {
  summaryDFfeature <- summaryDFfeature_other[summaryDFfeature_other$libName ==
                                               otherNamesPlot[x],]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = otherColours[x]) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = otherColours[x]) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_list_other[[x]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_list_other[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_list_other[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranFeat
ggObj2_combined_other <- lapply(seq_along(otherNames), function(x) {
  summaryDFranFeat <- summaryDFranFeat_other[summaryDFranFeat_other$libName ==
                                               otherNamesPlot[x],]
  ggplot(data = summaryDFranFeat,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranFeat,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = otherColours[x]) +
  geom_ribbon(data = summaryDFranFeat,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = otherColours[x]) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranFeat_list_other[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranFeat_list_other[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranFeat_list_other[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("ASY1 Cluster" * .(clusterLast) ~ "genes (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranFeat$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranLoc
ggObj3_combined_other <- lapply(seq_along(otherNames), function(x) {
  summaryDFranLoc <- summaryDFranLoc_other[summaryDFranLoc_other$libName ==
                                             otherNamesPlot[x],]
  ggplot(data = summaryDFranLoc,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranLoc,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = otherColours[x]) +
  geom_ribbon(data = summaryDFranLoc,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = otherColours[x]) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranLoc_list_other[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranLoc_list_other[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranLoc_list_other[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("Random loci (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranLoc$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## feature
ggObj1_combined_sRNA <- lapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[summaryDFfeature_sRNA$libName ==
                                               sRNANamesPlot[x],]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = sRNAColours[x]) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = sRNAColours[x]) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_list_sRNA[[x]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_list_sRNA[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_list_sRNA[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranFeat
ggObj2_combined_sRNA <- lapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFranFeat <- summaryDFranFeat_sRNA[summaryDFranFeat_sRNA$libName ==
                                               sRNANamesPlot[x],]
  ggplot(data = summaryDFranFeat,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranFeat,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = sRNAColours[x]) +
  geom_ribbon(data = summaryDFranFeat,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = sRNAColours[x]) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranFeat_list_sRNA[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranFeat_list_sRNA[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranFeat_list_sRNA[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("ASY1 Cluster" * .(clusterLast) ~ "genes (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranFeat$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranLoc
ggObj3_combined_sRNA <- lapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFranLoc <- summaryDFranLoc_sRNA[summaryDFranLoc_sRNA$libName ==
                                             sRNANamesPlot[x],]
  ggplot(data = summaryDFranLoc,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranLoc,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = sRNAColours[x]) +
  geom_ribbon(data = summaryDFranLoc,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = sRNAColours[x]) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranLoc_list_sRNA[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranLoc_list_sRNA[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranLoc_list_sRNA[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("Random loci (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranLoc$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## feature
ggObj1_combined_DNAmeth <- lapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[summaryDFfeature_DNAmeth$libName ==
                                                 DNAmethNamesPlot[x],]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = DNAmethColours[x]) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = DNAmethColours[x]) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_list_DNAmeth[[x]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_list_DNAmeth[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_list_DNAmeth[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranFeat
ggObj2_combined_DNAmeth <- lapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFranFeat <- summaryDFranFeat_DNAmeth[summaryDFranFeat_DNAmeth$libName ==
                                               DNAmethNamesPlot[x],]
  ggplot(data = summaryDFranFeat,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranFeat,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = DNAmethColours[x]) +
  geom_ribbon(data = summaryDFranFeat,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = DNAmethColours[x]) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranFeat_list_DNAmeth[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranFeat_list_DNAmeth[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranFeat_list_DNAmeth[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("ASY1 Cluster" * .(clusterLast) ~ "genes (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranFeat$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
## ranLoc
ggObj3_combined_DNAmeth <- lapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFranLoc <- summaryDFranLoc_DNAmeth[summaryDFranLoc_DNAmeth$libName ==
                                             DNAmethNamesPlot[x],]
  ggplot(data = summaryDFranLoc,
         mapping = aes(x = winNo,
                       y = mean,
                       group = libName),
        ) +
  geom_line(data = summaryDFranLoc,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = DNAmethColours[x]) +
  geom_ribbon(data = summaryDFranLoc,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = DNAmethColours[x]) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranLoc_list_DNAmeth[[x]])[1])-(downstream/binSize),
                              dim(summaryDFranLoc_list_DNAmeth[[x]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranLoc_list_DNAmeth[[x]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("Random loci (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranLoc$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
})
ggObjGA_combined <- grid.arrange(grobs = c(ggObj1_combined_log2ChIP,
                                           ggObj1_combined_other,
                                           ggObj1_combined_sRNA,
                                           ggObj1_combined_DNAmeth,
                                           ggObj2_combined_log2ChIP,
                                           ggObj2_combined_other,
                                           ggObj2_combined_sRNA,
                                           ggObj2_combined_DNAmeth,
                                           ggObj3_combined_log2ChIP,
                                           ggObj3_combined_other,
                                           ggObj3_combined_sRNA,
                                           ggObj3_combined_DNAmeth),
                                 layout_matrix = cbind(1:length(c(ChIPNames, otherNames, sRNANamesPlot, DNAmethNamesPlot)),
                                                       (length(c(ChIPNames, otherNames, sRNANamesPlot, DNAmethNamesPlot))+1):(length(c(ChIPNames, otherNames, sRNANamesPlot, DNAmethNamesPlot))*2),
                                                       ((length(c(ChIPNames, otherNames, sRNANamesPlot, DNAmethNamesPlot))*2)+1):(length(c(ChIPNames, otherNames, sRNANamesPlot, DNAmethNamesPlot))*3)))
ggsave(paste0(plotDir,
              "avgProfiles_around_",
              featureNamePlot, "_in_cluster", clusterNo,
              "_by_log2_ASY1_CS_Rep1_ChIP_control_in_", region,
              "_of_", featureName, "_v141119.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(ChIPNames, otherNames, sRNANamesPlot, DNAmethNamesPlot)), width = 21, limitsize = FALSE)
