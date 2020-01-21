#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage profiles with 95% CIs around CBF-encoding genes

# Usage:
# /applications/R/R-3.5.0/bin/Rscript CBF_genes_avgProfileRibbon.R ASY1_CS_Rep1_ChIP ASY1_CS both 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' 'CBF-encoding_genes' 3500 2000 2kb '2 kb' 20 20bp promoters

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#align <- "both"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#featureNamePlot <- "CBF-encoding_genes"
#bodyLength <- 3500
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "promoters"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
align <- args[3]
featureName <- unlist(strsplit(args[4],
                               split = ","))
featureNamePlot <- args[5]
bodyLength <- as.numeric(args[6])
upstream <- as.numeric(args[7])
downstream <- as.numeric(args[7])
flankName <- args[8]
flankNamePlot <- args[9]
binSize <- as.numeric(args[10])
binName <- args[11]
region <- args[12]

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

# Define plot titles
featureNamePlot <- featureNamePlot
ranFeatNamePlot <- "Random genes"
ranLocNamePlot <- "Random loci"

# Define feature start and end labels for plotting
if(grepl("genes", featureName)) {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Load CBFs
CBFs <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_manually_curated_gene_families/IWGSC_v1.1_cbf_representative_mRNA.gff3",
                   header = F)
 
# Load features
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
featureIDs <- sub(pattern = "\\.\\d+", replacement = "",
                  features$V9)

# Subset features to only those corresponding to CBFs
features_CBFs <- features[features$V1 %in% CBFs$V1 &
                          (features$V4 %in% CBFs$V4 |
                           features$V5 %in% CBFs$V5 ) &
                          features$V7 %in% CBFs$V7,]
# Get CBF IDs and their row indices in features
IDs <- sub(pattern = "\\.\\d+", replacement = "",
           features_CBFs$V9)
ID_indices <- which(featureIDs %in% IDs)
nonIDs <- featureIDs[!(featureIDs %in% IDs)]
# Function to randomly select feature IDs not present in IDs
ran_nonIDs_select <- function(nonIDsChr, n) {
  sample(x = nonIDsChr,
         size = n,
         replace = FALSE)
}
# Apply ran_nonIDs_select() function on a per-chromosome basis
# and create growing vector of feature IDs called ran_nonIDs
set.seed(9237452)
ran_nonIDs <- NULL
for(i in 1:length(levels(features$V1))) {
  IDsChr <- IDs[grepl(paste0("TraesCS", i), IDs)]
  nonIDsChr <- nonIDs[grepl(paste0("TraesCS", i), nonIDs)]
  ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                     n = length(IDsChr))
  ran_nonIDs <- c(ran_nonIDs, ran_nonIDsChr)
}
ran_nonID_indices <- which(featureIDs %in% ran_nonIDs)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]

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
log2ChIPNamesPlot <- c(
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
log2ChIPColours <- c(
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
      stop(paste0("controlNames[", x, "] is neither H3_input_SRR6350669 nor MNase_Rep1"))
    }
  }
})

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
                   "33-nt sRNAs",
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
sRNADirs <- sapply(seq_along(sRNANames), function(x) {
  if(sRNANames[x] %in% c("CS+_2_LIB18613_LDI16228")) {
    paste0("/home/ajt200/analysis/wheat/",
           sRNANamesDir[x], "/snakemake_sRNAseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("sRNANames[", x, "] is not compatible with the specified coverage matrix paths"))
  }
})

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
DNAmethDirs <- sapply(seq_along(DNAmethNames), function(x) {
  if(DNAmethNames[x] %in% c("BSseq_Rep8a_SRR6792678")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           DNAmethNamesDir[x],
           "/snakemake_BSseq/coverage/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("DNAmethNames[", x, "] is not compatible with the specified coverage matrix paths"))
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
log2ChIP_mats <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  list(
       # features 
       log2ChIP_featureMats[[x]][ID_indices,],
       # random features
       log2ChIP_featureMats[[x]][ran_nonID_indices,],
       # random loci
       log2ChIP_ranLocMats[[x]][ID_indices,]
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

names(summaryDFfeature_list_log2ChIP) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- summaryDFfeature_list_log2ChIP

# Define y-axis limits
ymin_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  min(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_lower))
})
ymax_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  max(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = log2ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = log2ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.9, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranFeat
ggObj2_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = log2ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = log2ChIPColours[x],
              alpha = 0.4) +
  scale_fill_manual(values = log2ChIPColours[x]) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranLoc
ggObj3_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = log2ChIPColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = log2ChIPColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_log2ChIP,
                                           ggObj2_combined_log2ChIP,
                                           ggObj3_combined_log2ChIP
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(log2ChIPNamesPlot)),
                                                       (length(c(log2ChIPNamesPlot))+1):(length(c(log2ChIPNamesPlot))*2),
                                                       ((length(c(log2ChIPNamesPlot))*2)+1):(length(c(log2ChIPNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_avgProfiles_around_",
              featureNamePlot, "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v200120.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(log2ChIPNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   ChIP_featureMats, ChIP_ranLocMats,
   control_featureMats, control_ranLocMats,
   log2ChIP_featureMats, log2ChIP_ranLocMats,
   log2ChIP_mats,
   wideDFfeature_list_log2ChIP,
   tidyDFfeature_list_log2ChIP,
   summaryDFfeature_list_log2ChIP,
   summaryDFfeature_log2ChIP
  ) 
gc()
#####



## other
# feature
other_featureMats <- mclapply(seq_along(otherNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    otherFile <- system(paste0("ls ", otherDirs[x],
                               otherNames[x],
                               "_MappedOn_wheat_v1.0*", align, "_sort_norm_",
                               featureName[y], "_matrix_bin", binName,
                               "_flank", flankName, ".tab"),
                        intern = T)
    as.matrix(read.table(otherFile,
                         header = F, skip = 3))
  })
}, mc.cores = length(otherNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
other_featureMats <- mclapply(seq_along(other_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, other_featureMats[[x]])
  } else {
    other_featureMats[[x]][[1]]
  }
}, mc.cores = length(other_featureMats))

# ranLoc
other_ranLocMats <- mclapply(seq_along(otherNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    otherFile <- system(paste0("ls ", otherDirs[x],
                               otherNames[x],
                               "_MappedOn_wheat_v1.0*", align, "_sort_norm_",
                               featureName[y], "_ranLoc_matrix_bin", binName,
                               "_flank", flankName, ".tab"),
                        intern = T)
    as.matrix(read.table(otherFile,
                         header = F, skip = 3))
  })
}, mc.cores = length(otherNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
other_ranLocMats <- mclapply(seq_along(other_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, other_ranLocMats[[x]])
  } else {
    other_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(other_ranLocMats))

# Add column names
for(x in seq_along(other_featureMats)) {
  colnames(other_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(other_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
other_mats <- mclapply(seq_along(other_featureMats), function(x) {
  list(
       # features 
       other_featureMats[[x]][ID_indices,],
       # random features
       other_featureMats[[x]][ran_nonID_indices,],
       # random loci
       other_ranLocMats[[x]][ID_indices,]
      ) 
}, mc.cores = length(other_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_other <- mclapply(seq_along(other_mats), function(x) {
  lapply(seq_along(other_mats[[x]]), function(y) {
      data.frame(window = colnames(other_mats[[x]][[y]]),
                 t(other_mats[[x]][[y]]))
  })
}, mc.cores = length(other_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_other  <- mclapply(seq_along(wideDFfeature_list_other), function(x) {
  lapply(seq_along(other_mats[[x]]), function(y) {
      gather(data  = wideDFfeature_list_other[[x]][[y]],
             key   = feature,
             value = coverage,
             -window)
  }) 
}, mc.cores = length(wideDFfeature_list_other))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_other)) {
  for(y in seq_along(other_mats[[x]])) {
      tidyDFfeature_list_other[[x]][[y]]$window <- factor(tidyDFfeature_list_other[[x]][[y]]$window,
                                                             levels = as.character(wideDFfeature_list_other[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_other  <- mclapply(seq_along(tidyDFfeature_list_other), function(x) {
  lapply(seq_along(other_mats[[x]]), function(y) {
      data.frame(window = as.character(wideDFfeature_list_other[[x]][[y]]$window),
                 n      = tapply(X     = tidyDFfeature_list_other[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_other[[x]][[y]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_other[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_other[[x]][[y]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_other[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_other[[x]][[y]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_other))

for(x in seq_along(summaryDFfeature_list_other)) {
  for(y in seq_along(other_mats[[x]])) {
      summaryDFfeature_list_other[[x]][[y]]$window <- factor(summaryDFfeature_list_other[[x]][[y]]$window,
                                                                levels = as.character(wideDFfeature_list_other[[x]][[y]]$window))
      summaryDFfeature_list_other[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_other[[x]][[y]])[1])
      summaryDFfeature_list_other[[x]][[y]]$sem <- summaryDFfeature_list_other[[x]][[y]]$sd/sqrt(summaryDFfeature_list_other[[x]][[y]]$n-1)
      summaryDFfeature_list_other[[x]][[y]]$CI_lower <- summaryDFfeature_list_other[[x]][[y]]$mean -
        qt(0.975, df = summaryDFfeature_list_other[[x]][[y]]$n-1)*summaryDFfeature_list_other[[x]][[y]]$sem
      summaryDFfeature_list_other[[x]][[y]]$CI_upper <- summaryDFfeature_list_other[[x]][[y]]$mean +
        qt(0.975, df = summaryDFfeature_list_other[[x]][[y]]$n-1)*summaryDFfeature_list_other[[x]][[y]]$sem
  }
}

names(summaryDFfeature_list_other) <- otherNamesPlot
summaryDFfeature_other <- summaryDFfeature_list_other

# Define y-axis limits
ymin_list_other <- lapply(seq_along(summaryDFfeature_other), function(x) {
  min(c(summaryDFfeature_other[[x]][[1]]$CI_lower,
        summaryDFfeature_other[[x]][[2]]$CI_lower,
        summaryDFfeature_other[[x]][[3]]$CI_lower))
})
ymax_list_other <- lapply(seq_along(summaryDFfeature_other), function(x) {
  max(c(summaryDFfeature_other[[x]][[1]]$CI_upper,
        summaryDFfeature_other[[x]][[2]]$CI_upper,
        summaryDFfeature_other[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = otherColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = otherColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[1]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 0.9, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranFeat
ggObj2_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = otherColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = otherColours[x],
              alpha = 0.4) +
  scale_fill_manual(values = otherColours[x]) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[2]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[2]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranLoc
ggObj3_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = otherColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = otherColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[3]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[3]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[3]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_other,
                                           ggObj2_combined_other,
                                           ggObj3_combined_other
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(otherNamesPlot)),
                                                       (length(c(otherNamesPlot))+1):(length(c(otherNamesPlot))*2),
                                                       ((length(c(otherNamesPlot))*2)+1):(length(c(otherNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "other_avgProfiles_around_",
              featureNamePlot, "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v200120.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(otherNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   other_featureMats, other_ranLocMats,
   other_mats,
   wideDFfeature_list_other,
   tidyDFfeature_list_other,
   summaryDFfeature_list_other,
   summaryDFfeature_other
  ) 
gc()
#####


## sRNA
# feature
sRNA_featureMats <- mclapply(seq_along(sRNAsizes), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(sRNADirs,
                                sRNANames,
                                "_MappedOn_wheat_v1.0_", align, "_", sRNAsizes[x], "_sort_norm_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(sRNAsizes))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
sRNA_featureMats <- mclapply(seq_along(sRNA_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, sRNA_featureMats[[x]])
  } else {
    sRNA_featureMats[[x]][[1]]
  }
}, mc.cores = length(sRNA_featureMats))

# ranLoc
sRNA_ranLocMats <- mclapply(seq_along(sRNAsizes), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(sRNADirs,
                                sRNANames,
                                "_MappedOn_wheat_v1.0_", align, "_", sRNAsizes[x], "_sort_norm_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(sRNAsizes))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
sRNA_ranLocMats <- mclapply(seq_along(sRNA_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, sRNA_ranLocMats[[x]])
  } else {
    sRNA_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(sRNA_ranLocMats))

# Add column names
for(x in seq_along(sRNA_featureMats)) {
  colnames(sRNA_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(sRNA_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
sRNA_mats <- mclapply(seq_along(sRNA_featureMats), function(x) {
  list(
       # features 
       sRNA_featureMats[[x]][ID_indices,],
       # random features
       sRNA_featureMats[[x]][ran_nonID_indices,],
       # random loci
       sRNA_ranLocMats[[x]][ID_indices,]
      ) 
}, mc.cores = length(sRNA_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_sRNA <- mclapply(seq_along(sRNA_mats), function(x) {
  lapply(seq_along(sRNA_mats[[x]]), function(y) {
      data.frame(window = colnames(sRNA_mats[[x]][[y]]),
                 t(sRNA_mats[[x]][[y]]))
  })
}, mc.cores = length(sRNA_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_sRNA  <- mclapply(seq_along(wideDFfeature_list_sRNA), function(x) {
  lapply(seq_along(sRNA_mats[[x]]), function(y) {
      gather(data  = wideDFfeature_list_sRNA[[x]][[y]],
             key   = feature,
             value = coverage,
             -window)
  }) 
}, mc.cores = length(wideDFfeature_list_sRNA))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_sRNA)) {
  for(y in seq_along(sRNA_mats[[x]])) {
      tidyDFfeature_list_sRNA[[x]][[y]]$window <- factor(tidyDFfeature_list_sRNA[[x]][[y]]$window,
                                                             levels = as.character(wideDFfeature_list_sRNA[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_sRNA  <- mclapply(seq_along(tidyDFfeature_list_sRNA), function(x) {
  lapply(seq_along(sRNA_mats[[x]]), function(y) {
      data.frame(window = as.character(wideDFfeature_list_sRNA[[x]][[y]]$window),
                 n      = tapply(X     = tidyDFfeature_list_sRNA[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_sRNA[[x]][[y]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_sRNA[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_sRNA[[x]][[y]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_sRNA[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_sRNA[[x]][[y]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_sRNA))

for(x in seq_along(summaryDFfeature_list_sRNA)) {
  for(y in seq_along(sRNA_mats[[x]])) {
      summaryDFfeature_list_sRNA[[x]][[y]]$window <- factor(summaryDFfeature_list_sRNA[[x]][[y]]$window,
                                                                levels = as.character(wideDFfeature_list_sRNA[[x]][[y]]$window))
      summaryDFfeature_list_sRNA[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_sRNA[[x]][[y]])[1])
      summaryDFfeature_list_sRNA[[x]][[y]]$sem <- summaryDFfeature_list_sRNA[[x]][[y]]$sd/sqrt(summaryDFfeature_list_sRNA[[x]][[y]]$n-1)
      summaryDFfeature_list_sRNA[[x]][[y]]$CI_lower <- summaryDFfeature_list_sRNA[[x]][[y]]$mean -
        qt(0.975, df = summaryDFfeature_list_sRNA[[x]][[y]]$n-1)*summaryDFfeature_list_sRNA[[x]][[y]]$sem
      summaryDFfeature_list_sRNA[[x]][[y]]$CI_upper <- summaryDFfeature_list_sRNA[[x]][[y]]$mean +
        qt(0.975, df = summaryDFfeature_list_sRNA[[x]][[y]]$n-1)*summaryDFfeature_list_sRNA[[x]][[y]]$sem
  }
}

names(summaryDFfeature_list_sRNA) <- sRNANamesPlot
summaryDFfeature_sRNA <- summaryDFfeature_list_sRNA

# Define y-axis limits
ymin_list_sRNA <- lapply(seq_along(summaryDFfeature_sRNA), function(x) {
  min(c(summaryDFfeature_sRNA[[x]][[1]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[2]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[3]]$CI_lower))
})
ymax_list_sRNA <- lapply(seq_along(summaryDFfeature_sRNA), function(x) {
  max(c(summaryDFfeature_sRNA[[x]][[1]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[2]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = sRNAColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = sRNAColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[1]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 0.9, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranFeat
ggObj2_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = sRNAColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = sRNAColours[x],
              alpha = 0.4) +
  scale_fill_manual(values = sRNAColours[x]) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[2]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[2]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranLoc
ggObj3_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = sRNAColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = sRNAColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[3]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[3]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[3]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_sRNA,
                                           ggObj2_combined_sRNA,
                                           ggObj3_combined_sRNA
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(sRNANamesPlot)),
                                                       (length(c(sRNANamesPlot))+1):(length(c(sRNANamesPlot))*2),
                                                       ((length(c(sRNANamesPlot))*2)+1):(length(c(sRNANamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "sRNA_avgProfiles_around_",
              featureNamePlot, "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v200120.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(sRNANamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   sRNA_featureMats, sRNA_ranLocMats,
   sRNA_mats,
   wideDFfeature_list_sRNA,
   tidyDFfeature_list_sRNA,
   summaryDFfeature_list_sRNA,
   summaryDFfeature_sRNA
  ) 
gc()
#####


## DNAmeth
# feature
DNAmeth_featureMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs,
                                DNAmethNames,
                                "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_", DNAmethContexts[x], "_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
DNAmeth_featureMats <- mclapply(seq_along(DNAmeth_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, DNAmeth_featureMats[[x]])
  } else {
    DNAmeth_featureMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_featureMats))

# ranLoc
DNAmeth_ranLocMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs,
                                DNAmethNames,
                                "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_", DNAmethContexts[x], "_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
DNAmeth_ranLocMats <- mclapply(seq_along(DNAmeth_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, DNAmeth_ranLocMats[[x]])
  } else {
    DNAmeth_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_ranLocMats))

# Add column names
for(x in seq_along(DNAmeth_featureMats)) {
  colnames(DNAmeth_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(DNAmeth_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
DNAmeth_mats <- mclapply(seq_along(DNAmeth_featureMats), function(x) {
  list(
       # features 
       DNAmeth_featureMats[[x]][ID_indices,],
       # random features
       DNAmeth_featureMats[[x]][ran_nonID_indices,],
       # random loci
       DNAmeth_ranLocMats[[x]][ID_indices,]
      ) 
}, mc.cores = length(DNAmeth_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_DNAmeth <- mclapply(seq_along(DNAmeth_mats), function(x) {
  lapply(seq_along(DNAmeth_mats[[x]]), function(y) {
      data.frame(window = colnames(DNAmeth_mats[[x]][[y]]),
                 t(DNAmeth_mats[[x]][[y]]))
  })
}, mc.cores = length(DNAmeth_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_DNAmeth  <- mclapply(seq_along(wideDFfeature_list_DNAmeth), function(x) {
  lapply(seq_along(DNAmeth_mats[[x]]), function(y) {
      gather(data  = wideDFfeature_list_DNAmeth[[x]][[y]],
             key   = feature,
             value = coverage,
             -window)
  }) 
}, mc.cores = length(wideDFfeature_list_DNAmeth))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_DNAmeth)) {
  for(y in seq_along(DNAmeth_mats[[x]])) {
      tidyDFfeature_list_DNAmeth[[x]][[y]]$window <- factor(tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                                                             levels = as.character(wideDFfeature_list_DNAmeth[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_DNAmeth  <- mclapply(seq_along(tidyDFfeature_list_DNAmeth), function(x) {
  lapply(seq_along(DNAmeth_mats[[x]]), function(y) {
      data.frame(window = as.character(wideDFfeature_list_DNAmeth[[x]][[y]]$window),
                 n      = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_DNAmeth))

for(x in seq_along(summaryDFfeature_list_DNAmeth)) {
  for(y in seq_along(DNAmeth_mats[[x]])) {
      summaryDFfeature_list_DNAmeth[[x]][[y]]$window <- factor(summaryDFfeature_list_DNAmeth[[x]][[y]]$window,
                                                                levels = as.character(wideDFfeature_list_DNAmeth[[x]][[y]]$window))
      summaryDFfeature_list_DNAmeth[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_DNAmeth[[x]][[y]])[1])
      summaryDFfeature_list_DNAmeth[[x]][[y]]$sem <- summaryDFfeature_list_DNAmeth[[x]][[y]]$sd/sqrt(summaryDFfeature_list_DNAmeth[[x]][[y]]$n-1)
      summaryDFfeature_list_DNAmeth[[x]][[y]]$CI_lower <- summaryDFfeature_list_DNAmeth[[x]][[y]]$mean -
        qt(0.975, df = summaryDFfeature_list_DNAmeth[[x]][[y]]$n-1)*summaryDFfeature_list_DNAmeth[[x]][[y]]$sem
      summaryDFfeature_list_DNAmeth[[x]][[y]]$CI_upper <- summaryDFfeature_list_DNAmeth[[x]][[y]]$mean +
        qt(0.975, df = summaryDFfeature_list_DNAmeth[[x]][[y]]$n-1)*summaryDFfeature_list_DNAmeth[[x]][[y]]$sem
  }
}

names(summaryDFfeature_list_DNAmeth) <- DNAmethNamesPlot
summaryDFfeature_DNAmeth <- summaryDFfeature_list_DNAmeth

# Define y-axis limits
ymin_list_DNAmeth <- lapply(seq_along(summaryDFfeature_DNAmeth), function(x) {
  min(c(summaryDFfeature_DNAmeth[[x]][[1]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[2]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[3]]$CI_lower))
})
ymax_list_DNAmeth <- lapply(seq_along(summaryDFfeature_DNAmeth), function(x) {
  max(c(summaryDFfeature_DNAmeth[[x]][[1]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[2]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = DNAmethColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = DNAmethColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[1]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 0.9, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranFeat
ggObj2_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = DNAmethColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = DNAmethColours[x],
              alpha = 0.4) +
  scale_fill_manual(values = DNAmethColours[x]) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[2]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[2]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranLoc
ggObj3_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = DNAmethColours[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = DNAmethColours[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[3]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[3]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[3]])[1])-(downstream/binSize)),
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
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_DNAmeth,
                                           ggObj2_combined_DNAmeth,
                                           ggObj3_combined_DNAmeth
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(DNAmethNamesPlot)),
                                                       (length(c(DNAmethNamesPlot))+1):(length(c(DNAmethNamesPlot))*2),
                                                       ((length(c(DNAmethNamesPlot))*2)+1):(length(c(DNAmethNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "DNAmeth_avgProfiles_around_",
              featureNamePlot, "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v200120.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(DNAmethNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   DNAmeth_featureMats, DNAmeth_ranLocMats,
   DNAmeth_mats,
   wideDFfeature_list_DNAmeth,
   tidyDFfeature_list_DNAmeth,
   summaryDFfeature_list_DNAmeth,
   summaryDFfeature_DNAmeth
  ) 
gc()
#####



# varietal SNPclasses
SNPclassNames <- c(
                   "all",
                   "missense_variant",
                   "synonymous_variant",
                   "HIGH",
                   "MODERATE",
                   "LOW",
                   "MODIFIER",
                   "upstream_gene_variant",
                   "downstream_gene_variant",
                   "intron_variant",
                   "intergenic",
                   "transition",
                   "transversion"
                  )
SNPclassNamesPlot <- c(
                       "All SNPs",
                       "Missense SNPs",
                       "Synonymous SNPs",
                       "High-impact SNPs",
                       "Moderate-impact SNPs",
                       "Low-impact SNPs",
                       "Modifier-impact SNPs",
                       "SNPs upstream of a gene",
                       "SNPs downstream of a gene",
                       "Intronic SNPs",
                       "Intergenic SNPs",
                       "Transitions",
                       "Transversions"
                      )

# feature
SNPclass_featureMats <- mclapply(seq_along(SNPclassNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                SNPclassNames[x],
                                "_SNPs_around_", featureName[y],
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(SNPclassNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
SNPclass_featureMats <- mclapply(seq_along(SNPclass_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, SNPclass_featureMats[[x]])
  } else {
    SNPclass_featureMats[[x]][[1]]
  }
}, mc.cores = length(SNPclass_featureMats))

# ranLoc
SNPclass_ranLocMats <- mclapply(seq_along(SNPclassNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                SNPclassNames[x],
                                "_SNPs_around_", featureName[y],
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(SNPclassNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
SNPclass_ranLocMats <- mclapply(seq_along(SNPclass_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, SNPclass_ranLocMats[[x]])
  } else {
    SNPclass_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(SNPclass_ranLocMats))

# Add column names
for(x in seq_along(SNPclass_featureMats)) {
  colnames(SNPclass_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(SNPclass_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
SNPclass_mats <- mclapply(seq_along(SNPclass_featureMats), function(x) {
  list(
       # features 
       SNPclass_featureMats[[x]][ID_indices,],
       # random features
       SNPclass_featureMats[[x]][ran_nonID_indices,],
       # random loci
       SNPclass_ranLocMats[[x]][ID_indices,]
      ) 
}, mc.cores = length(SNPclass_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_SNPclass <- mclapply(seq_along(SNPclass_mats), function(x) {
  lapply(seq_along(SNPclass_mats[[x]]), function(y) {
      data.frame(window = colnames(SNPclass_mats[[x]][[y]]),
                 t(SNPclass_mats[[x]][[y]]))
  })
}, mc.cores = length(SNPclass_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_SNPclass  <- mclapply(seq_along(wideDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats[[x]]), function(y) {
      gather(data  = wideDFfeature_list_SNPclass[[x]][[y]],
             key   = feature,
             value = coverage,
             -window)
  }) 
}, mc.cores = length(wideDFfeature_list_SNPclass))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_SNPclass)) {
  for(y in seq_along(SNPclass_mats[[x]])) {
      tidyDFfeature_list_SNPclass[[x]][[y]]$window <- factor(tidyDFfeature_list_SNPclass[[x]][[y]]$window,
                                                             levels = as.character(wideDFfeature_list_SNPclass[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_SNPclass  <- mclapply(seq_along(tidyDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats[[x]]), function(y) {
      data.frame(window = as.character(wideDFfeature_list_SNPclass[[x]][[y]]$window),
                 n      = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_SNPclass))

for(x in seq_along(summaryDFfeature_list_SNPclass)) {
  for(y in seq_along(SNPclass_mats[[x]])) {
      summaryDFfeature_list_SNPclass[[x]][[y]]$window <- factor(summaryDFfeature_list_SNPclass[[x]][[y]]$window,
                                                                levels = as.character(wideDFfeature_list_SNPclass[[x]][[y]]$window))
      summaryDFfeature_list_SNPclass[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_SNPclass[[x]][[y]])[1])
      summaryDFfeature_list_SNPclass[[x]][[y]]$sem <- summaryDFfeature_list_SNPclass[[x]][[y]]$sd/sqrt(summaryDFfeature_list_SNPclass[[x]][[y]]$n-1)
      summaryDFfeature_list_SNPclass[[x]][[y]]$CI_lower <- summaryDFfeature_list_SNPclass[[x]][[y]]$mean -
        qt(0.975, df = summaryDFfeature_list_SNPclass[[x]][[y]]$n-1)*summaryDFfeature_list_SNPclass[[x]][[y]]$sem
      summaryDFfeature_list_SNPclass[[x]][[y]]$CI_upper <- summaryDFfeature_list_SNPclass[[x]][[y]]$mean +
        qt(0.975, df = summaryDFfeature_list_SNPclass[[x]][[y]]$n-1)*summaryDFfeature_list_SNPclass[[x]][[y]]$sem
  }
}

names(summaryDFfeature_list_SNPclass) <- SNPclassNamesPlot
summaryDFfeature_SNPclass <- summaryDFfeature_list_SNPclass

# Define y-axis limits
ymin_list_SNPclass <- lapply(seq_along(summaryDFfeature_SNPclass), function(x) {
  min(c(summaryDFfeature_SNPclass[[x]][[1]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[2]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[3]]$CI_lower))
})
ymax_list_SNPclass <- lapply(seq_along(summaryDFfeature_SNPclass), function(x) {
  max(c(summaryDFfeature_SNPclass[[x]][[1]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[2]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = c(log2ChIPColours, otherColours)[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = c(log2ChIPColours, otherColours)[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = c(log2ChIPColours, otherColours)[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.9, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranFeat
ggObj2_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = c(log2ChIPColours, otherColours)[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = c(log2ChIPColours, otherColours)[x],
              alpha = 0.4) +
  scale_fill_manual(values = c(log2ChIPColours, otherColours)[x]) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[2]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[2]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = c(log2ChIPColours, otherColours)[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranLoc
ggObj3_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = c(log2ChIPColours, otherColours)[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = c(log2ChIPColours, otherColours)[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[3]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[3]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[3]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = c(log2ChIPColours, otherColours)[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_SNPclass,
                                           ggObj2_combined_SNPclass,
                                           ggObj3_combined_SNPclass
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(SNPclassNamesPlot)),
                                                       (length(c(SNPclassNamesPlot))+1):(length(c(SNPclassNamesPlot))*2),
                                                       ((length(c(SNPclassNamesPlot))*2)+1):(length(c(SNPclassNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "SNPclass_avgProfiles_around_",
              featureNamePlot, "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v200120.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(SNPclassNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   SNPclass_featureMats, SNPclass_ranLocMats,
   SNPclass_mats,
   wideDFfeature_list_SNPclass,
   tidyDFfeature_list_SNPclass,
   summaryDFfeature_list_SNPclass,
   summaryDFfeature_SNPclass
  ) 
gc()
#####



# TE superfams
superfamCodes <- c("RLG",
                   "RLC",
                   "RLX",
                   "RIX",
                   "SIX",
                   "DTC",
                   "DTM",
                   "DTX",
                   "DTH",
                   "DMI",
                   "DTT",
                   "DXX",
                   "DTA",
                   "DHH",
                   "XXX")
superfamNames <- c("Gypsy_LTR",
                   "Copia_LTR",
                   "Unclassified_LTR",
                   "LINE",
                   "SINE",
                   "CACTA",
                   "Mutator",
                   "Unclassified_with_TIRs",
                   "Harbinger",
                   "MITE",
                   "Mariner",
                   "Unclassified_class_2",
                   "hAT",
                   "Helitrons",
                   "Unclassified_repeats")
superfamNamesPlot <- c("Gypsy LTR",
                       "Copia LTR",
                       "Unclassified LTR",
                       "LINE",
                       "SINE",
                       "CACTA",
                       "Mutator",
                       "Unclassified with TIRs",
                       "Harbinger",
                       "MITE",
                       "Mariner",
                       "Unclassified class 2",
                       "hAT",
                       "Helitrons",
                       "Unclassified repeats")

# feature
superfam_featureMats <- mclapply(seq_along(superfamNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                superfamNames[x], "_", superfamCodes[x],
                                "_around_", featureName[y],
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(superfamNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
superfam_featureMats <- mclapply(seq_along(superfam_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, superfam_featureMats[[x]])
  } else {
    superfam_featureMats[[x]][[1]]
  }
}, mc.cores = length(superfam_featureMats))

# ranLoc
superfam_ranLocMats <- mclapply(seq_along(superfamNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                superfamNames[x], "_", superfamCodes[x],
                                "_around_", featureName[y],
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(superfamNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
superfam_ranLocMats <- mclapply(seq_along(superfam_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, superfam_ranLocMats[[x]])
  } else {
    superfam_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(superfam_ranLocMats))

# Add column names
for(x in seq_along(superfam_featureMats)) {
  colnames(superfam_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(superfam_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
superfam_mats <- mclapply(seq_along(superfam_featureMats), function(x) {
  list(
       # features 
       superfam_featureMats[[x]][ID_indices,],
       # random features
       superfam_featureMats[[x]][ran_nonID_indices,],
       # random loci
       superfam_ranLocMats[[x]][ID_indices,]
      ) 
}, mc.cores = length(superfam_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_superfam <- mclapply(seq_along(superfam_mats), function(x) {
  lapply(seq_along(superfam_mats[[x]]), function(y) {
      data.frame(window = colnames(superfam_mats[[x]][[y]]),
                 t(superfam_mats[[x]][[y]]))
  })
}, mc.cores = length(superfam_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_superfam  <- mclapply(seq_along(wideDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats[[x]]), function(y) {
      gather(data  = wideDFfeature_list_superfam[[x]][[y]],
             key   = feature,
             value = coverage,
             -window)
  }) 
}, mc.cores = length(wideDFfeature_list_superfam))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_superfam)) {
  for(y in seq_along(superfam_mats[[x]])) {
      tidyDFfeature_list_superfam[[x]][[y]]$window <- factor(tidyDFfeature_list_superfam[[x]][[y]]$window,
                                                             levels = as.character(wideDFfeature_list_superfam[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_superfam  <- mclapply(seq_along(tidyDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats[[x]]), function(y) {
      data.frame(window = as.character(wideDFfeature_list_superfam[[x]][[y]]$window),
                 n      = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_superfam))

for(x in seq_along(summaryDFfeature_list_superfam)) {
  for(y in seq_along(superfam_mats[[x]])) {
      summaryDFfeature_list_superfam[[x]][[y]]$window <- factor(summaryDFfeature_list_superfam[[x]][[y]]$window,
                                                                levels = as.character(wideDFfeature_list_superfam[[x]][[y]]$window))
      summaryDFfeature_list_superfam[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_superfam[[x]][[y]])[1])
      summaryDFfeature_list_superfam[[x]][[y]]$sem <- summaryDFfeature_list_superfam[[x]][[y]]$sd/sqrt(summaryDFfeature_list_superfam[[x]][[y]]$n-1)
      summaryDFfeature_list_superfam[[x]][[y]]$CI_lower <- summaryDFfeature_list_superfam[[x]][[y]]$mean -
        qt(0.975, df = summaryDFfeature_list_superfam[[x]][[y]]$n-1)*summaryDFfeature_list_superfam[[x]][[y]]$sem
      summaryDFfeature_list_superfam[[x]][[y]]$CI_upper <- summaryDFfeature_list_superfam[[x]][[y]]$mean +
        qt(0.975, df = summaryDFfeature_list_superfam[[x]][[y]]$n-1)*summaryDFfeature_list_superfam[[x]][[y]]$sem
  }
}

names(summaryDFfeature_list_superfam) <- superfamNamesPlot
summaryDFfeature_superfam <- summaryDFfeature_list_superfam

# Define y-axis limits
ymin_list_superfam <- lapply(seq_along(summaryDFfeature_superfam), function(x) {
  min(c(summaryDFfeature_superfam[[x]][[1]]$CI_lower,
        summaryDFfeature_superfam[[x]][[2]]$CI_lower,
        summaryDFfeature_superfam[[x]][[3]]$CI_lower))
})
ymax_list_superfam <- lapply(seq_along(summaryDFfeature_superfam), function(x) {
  max(c(summaryDFfeature_superfam[[x]][[1]]$CI_upper,
        summaryDFfeature_superfam[[x]][[2]]$CI_upper,
        summaryDFfeature_superfam[[x]][[3]]$CI_upper))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = c(log2ChIPColours, sRNAColours)[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = c(log2ChIPColours, sRNAColours)[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = c(log2ChIPColours, sRNAColours)[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.9, size = 30)) +
  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranFeat
ggObj2_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = c(log2ChIPColours, sRNAColours)[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = c(log2ChIPColours, sRNAColours)[x],
              alpha = 0.4) +
  scale_fill_manual(values = c(log2ChIPColours, sRNAColours)[x]) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[2]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[2]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[2]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = c(log2ChIPColours, sRNAColours)[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranLoc
ggObj3_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = 1)
        ) +
  geom_line(data = summaryDFfeature,
            colour = c(log2ChIPColours, sRNAColours)[x],
            size = 1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper),
              fill = c(log2ChIPColours, sRNAColours)[x],
              alpha = 0.4) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[3]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[3]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[3]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = c(log2ChIPColours, sRNAColours)[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_superfam,
                                           ggObj2_combined_superfam,
                                           ggObj3_combined_superfam
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(superfamNamesPlot)),
                                                       (length(c(superfamNamesPlot))+1):(length(c(superfamNamesPlot))*2),
                                                       ((length(c(superfamNamesPlot))*2)+1):(length(c(superfamNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "superfam_avgProfiles_around_",
              featureNamePlot, "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v200120.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(superfamNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   superfam_featureMats, superfam_ranLocMats,
   superfam_mats,
   wideDFfeature_list_superfam,
   tidyDFfeature_list_superfam,
   summaryDFfeature_list_superfam,
   summaryDFfeature_superfam
  ) 
gc()
#####

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_log2ChIP,
                                           ggObj1_combined_other,
                                           ggObj1_combined_sRNA,
                                           ggObj1_combined_DNAmeth,
                                           ggObj1_combined_SNPclass,
                                           ggObj1_combined_superfam,
                                           ggObj2_combined_log2ChIP,
                                           ggObj2_combined_other,
                                           ggObj2_combined_sRNA,
                                           ggObj2_combined_DNAmeth,
                                           ggObj2_combined_SNPclass,
                                           ggObj2_combined_superfam,
                                           ggObj3_combined_log2ChIP,
                                           ggObj3_combined_other,
                                           ggObj3_combined_sRNA,
                                           ggObj3_combined_DNAmeth,
                                           ggObj3_combined_SNPclass,
                                           ggObj3_combined_superfam
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot)),
                                                       (length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*2),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*2)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "combined_avgProfiles_around_",
              featureNamePlot, "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v200120.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot)), width = 21, limitsize = FALSE)

