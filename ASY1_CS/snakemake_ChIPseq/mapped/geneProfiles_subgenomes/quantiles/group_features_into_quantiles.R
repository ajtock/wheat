#!/applications/R/R-3.5.0/bin/Rscript

#
# Divide features into quantiles based on mean log2(libName ChIP/control)
# in a given feature region (e.g., promoters).
# Extract and save feature IDs for each quantile for further analyses
# (e.g., GO enrichment and average + 95% CI profile plotting).
# Calculate mean winName-scaled recombination rate (cM/Mb) from
# the promoter to the terminator of each feature.
# Plot feature quantile heatmaps for various genomics datasets
# Plot feature quantile recombination rate densities in a
# heat map or violin plot
#

# Usage:
# /applications/R/R-3.4.0/bin/Rscript group_features_into_quantiles.R ASY1_CS_Rep1_ChIP ASY1_CS both genes_in_Agenome_genomewide 3500 2000 2kb '2 kb' 20 20bp promoters 4 100kb 1

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

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
align <- args[3]
featureName <- args[4]
bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
flankNamePlot <- args[8]
binSize <- as.numeric(args[9])
binName <- args[10]
region <- args[11]
quantiles <- as.numeric(args[12])
winName <- args[13]
minMarkerDist <-  as.numeric(args[14])

library(GenomicRanges)
library(dplyr)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

outDir <- paste0("quantiles_by_log2_", libName,
                 "_control_in_", region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load ChIP matrix
mat1 <- lapply(seq_along(featureName), function(y) {
  as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/",
                              dirName,
                              "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                              libName,
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName[y], "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
})
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
if(length(featureName) == 3) {
 mat1 <- do.call(rbind, mat1)
} else {
 mat1 <- mat1[[1]]
}

# Load control matrices
controlNames <- c(
                  "H3_input_SRR6350669",
                  "MNase_Rep1"
                 )
controlNamesDir <- c(
                     "input",
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
controlmats <- mclapply(seq_along(controlNames), function(x) {
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
controlmats <- mclapply(seq_along(controlmats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, controlmats[[x]])
  } else {
    controlmats[[x]][[1]]
  }
}, mc.cores = length(controlmats))

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for ChIP matrix depending on library
log2ChIPmat <- if(libName %in% c(
                                 "ASY1_CS_Rep1_ChIP",
                                 "DMC1_Rep1_ChIP",
                                 "H3K4me3_ChIP_SRR6350668",
                                 "H3K27me3_ChIP_SRR6350666",
                                 "H3K36me3_ChIP_SRR6350670",
                                 "H3K9ac_ChIP_SRR6350667",
                                 #"DNaseI_Rep1_SRR8447247",
                                 "H3K4me1_Rep1_ChIP_SRR8126618",
                                 "H3K27ac_Rep1_ChIP_SRR8126621"
                                )) {
  print(paste0(libName, " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
  log2((mat1+1)/(controlmats[[1]]+1))
} else {
  print(paste0(libName, " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
  log2((mat1+1)/(controlmats[[2]]+1))
}

# Extract region for ordering of features (adjust promoter/terminator size as necessary)
if( region == "promoters" ) {
  log2ChIPmatRegion <- log2ChIPmat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  log2ChIPmatRegion <- log2ChIPmat[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
} else if ( region == "bodies" ) {
  log2ChIPmatRegion <- log2ChIPmat[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else {
  print("The region name provided does not match 'promoters', 'terminators', or 'bodies'")
}
log2ChIPmatRegionRowMeans <- rowMeans(log2ChIPmatRegion, na.rm = T)
log2ChIPmatRegionRowMeansSorted <- sort.int(log2ChIPmatRegionRowMeans,
                                            decreasing = T,
                                            index.return = T,
                                            na.last = T)
log2ChIPmatRegionSorted <- log2ChIPmatRegion[sort.int(log2ChIPmatRegionRowMeans,
                                                      decreasing = T,
                                                      index.return = T,
                                                      na.last = T)$ix,]
log2ChIPmatSorted <- log2ChIPmat[sort.int(log2ChIPmatRegionRowMeans,
                                          decreasing = T,
                                          index.return = T,
                                          na.last = T)$ix,]
## Replace NAs in log2ChIPmatRegion with 0
#log2ChIPmatRegion[which(is.na(log2ChIPmatRegion))] <- 0

# Load features 
features <- lapply(seq_along(featureName), function(y) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                    substring(featureName[y], first = 10), ".gff3"),
             colClasses = c(NA,
                            rep("NULL", 2),
                            rep(NA, 2),
                            "NULL", NA, "NULL", NA),
             header = F)
})
if(length(featureName) == 3) {
  features <- do.call(rbind, features)
} else {
  features <- features[[1]]
}
colnames(features) <- c("chr", "start", "end", "strand", "featureID")
featuresGR <- GRanges(seqnames = features$chr,
                      ranges = IRanges(start = features$start,
                                       end = features$end),
                      strand = features$strand,
                      featureID = features$featureID)
# Extend feature boundaries to include promoters and terminators for calculation of
# winName-scaled recombination rate
featuresGR_ext <- GRanges(seqnames = seqnames(featuresGR),
                          ranges = IRanges(start = start(featuresGR)-1000,
                                           end = end(featuresGR)+1000),
                          strand = strand(featuresGR),
                          featureID = featuresGR$featureID)
print(featuresGR_ext)

# Convert windowed recombination rate into GRanges
cMMb <- read.table(paste0(
                   "/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/",
                   "cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                   as.character(minMarkerDist), "bp_", winName, ".txt"))
cMMbGR <- GRanges(seqnames = cMMb$chr,
                  ranges = IRanges(start = cMMb$windowStart,
                                   end = cMMb$windowEnd),
                  strand = "*",
                  cMMb = cMMb$cMMb)

# Obtain winName-scaled cMMb values for each feature between promoter and terminator
# Where features overlap more than one winName window, calculate mean cMMb
feature_cMMb_overlaps <- findOverlaps(query = featuresGR_ext,
                                      subject = cMMbGR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
feature_cMMb_overlapsList <- lapply(seq_along(featuresGR_ext), function(x) {
  subjectHits(feature_cMMb_overlaps)[queryHits(feature_cMMb_overlaps) == x]
})
## OR
#feature_cMMb_overlapsList <- getOverlaps(coordinates = featuresGR_ext,
#                                         segments = cMMbGR,
#                                         overlapType = "overlapping",
#                                         whichOverlaps = TRUE,
#                                         ignoreStrand = TRUE)
feature_cMMb <- sapply(feature_cMMb_overlapsList,
                       function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
featuresGR <- GRanges(featuresGR,
                      featureID = featuresGR$featureID,
                      log2ChIPmatRegionRowMeans = log2ChIPmatRegionRowMeans,
                      cMMb = feature_cMMb)

# Divide features into quantiles based on decreasing log2ChIPmatRegionRowMeans
featuresDF <- data.frame(featuresGR,
                         quantile = as.character(""),
                         stringsAsFactors = F)
featuresDF$log2ChIPmatRegionRowMeans[which(is.na(featuresDF$log2ChIPmatRegionRowMeans))] <- 0
for(j in 1:quantiles) {
  if(j < quantiles) {
    featuresDF[ percent_rank(featuresDF$log2ChIPmatRegionRowMeans) <= (quantiles/quantiles)-((j-1)/quantiles) &
                percent_rank(featuresDF$log2ChIPmatRegionRowMeans) >  (quantiles/quantiles)-(j/quantiles), ]$quantile <- paste0("Quantile ", j)
  } else {
    featuresDF[ percent_rank(featuresDF$log2ChIPmatRegionRowMeans) <= (quantiles/quantiles)-((j-1)/quantiles) &
                percent_rank(featuresDF$log2ChIPmatRegionRowMeans) >= (quantiles/quantiles)-(j/quantiles), ]$quantile <- paste0("Quantile ", j)
  }
  write.table(featuresDF[featuresDF$quantile == paste0("Quantile ", j),],
              file = paste0(outDir,
                            "quantile", j, "_of_", quantiles,
                            "_by_log2_", libName, "_control_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}



  }
  if(j > 1 & j < quantiles) {
    quantilejFeatures <- featuresDF_ordered[(quantilesCum[j]):(quantilesCum[j+1]-1),]
    print(paste0("condition 2: quantile ", j))
    print(dim(quantilejFeatures))
  }
  if(j == quantiles) {
    quantilejFeatures <- featuresDF_ordered[(quantilesCum[j]):(quantilesCum[length(quantilesCum)]),]
    print(paste0("condition 3: quantile ", j))
    print(dim(quantilejFeatures))
  }
  quantilejFeatures <- data.frame(quantilejFeatures,
                                  quantile = paste0("Quantile ", j))
  if(j <= quantiles) {
    stats <- data.frame(quantile = as.integer(j),
                        n = as.integer(dim(quantilejFeatures)[1]),
                        mean_width = as.integer(round(mean(quantilejFeatures$width, na.rm = T))),
                        total_width = as.integer(sum(quantilejFeatures$width, na.rm = T)),
                        mean_log2ChIPmatRegionRowMeans = as.numeric(mean(quantilejFeatures$log2ChIPmatRegionRowMeans, na.rm = T)))
    quantilesStats <- rbind(quantilesStats, stats)
    write.table(quantilejFeatures,
                file = paste0(outDir,
                              "quantile", j, "_of_", quantiles,
                              "_by_log2_", libName, "_control_in_",
                              region, "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
}
write.table(quantilesStats,
            file = paste0(outDir,
                          "summary_", quantiles, "quantiles_by_log2_", libName, "_control_in_",
                          region, "_of_",
                          substring(featureName[1][1], first = 1, last = 5), "_in_",
                          paste0(substring(featureName, first = 10, last = 16),
                                 collapse = "_"), "_",
                          substring(featureName[1][1], first = 18), ".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)


  write.table(featureIDsClusterList[[k]],
              file = paste0(outDir, "quantile", as.character(k), "_of_", as.character(quantiles),
                            "_by_log2_", libName, "_control_in_",
                            region, "_of_", featureName, ".txt"),
              quote = F, row.names = F, col.names = F)
# Order features in each quantile by decreasing log2ChIPmatRegion levels
# to define "row_order" for heatmaps
combineRowOrders <- function(quantile_bool_list) {
  do.call("c", lapply(quantile_bool_list, function(x) {
    quantile_log2ChIPmatRegionRowMeans <- rowMeans(log2ChIPmatRegion[x,], na.rm = T)
    which(x)[order(quantile_log2ChIPmatRegionRowMeans, decreasing = T)]
  }))
}
row_order <- combineRowOrders(quantile_bool_list =
  lapply(seq_along(1:quantiles), function(k) { 
    km$quantile == paste0("Cluster ", k)
  })
)
# Order feature IDs in each quantile by decreasing log2ChIPmatRegion levels
# for use in GO term enrichment analysis
listCombineRowOrders <- function(quantile_bool_list) {
  do.call(list, lapply(quantile_bool_list, function(x) {
    quantile_log2ChIPmatRegionRowMeans <- rowMeans(log2ChIPmatRegion[x,], na.rm = T)
    which(x)[order(quantile_log2ChIPmatRegionRowMeans, decreasing = T)]
  }))
}
featureIndicesList <- listCombineRowOrders(quantile_bool_list =
  lapply(seq_along(1:quantiles), function(k) {
    km$quantile == paste0("Cluster ", k)
  })
)
# Alternatively, with original ordering:
## Get feature indices for each quantile
#featureIndicesList <- lapply(seq_along(1:quantiles), function(k) {
#  which(km$quantile == paste0("Cluster ", k))
#})

# Load features 
features <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                              substring(featureName, first = 10), ".gff3"),
                       header = F)
# Separate into quantiles
featuresClusterList <- lapply(seq_along(1:quantiles), function(k) {
  features[featureIndicesList[[k]],]
})
featureIDsClusterList <- lapply(seq_along(1:quantiles), function(k) {
  sub(pattern = "\\.\\d+", replacement = "",
      x = as.vector(features[featureIndicesList[[k]],]$V9))
})
sapply(seq_along(featureIDsClusterList), function(k) {
  write.table(featureIDsClusterList[[k]],
              file = paste0(outDir, "quantile", as.character(k), "_of_", as.character(quantiles),
                            "_by_log2_", libName, "_control_in_",
                            region, "_of_", featureName, ".txt"),
              quote = F, row.names = F, col.names = F)
})

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
ChIPmats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x],
                              ChIPNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))
othermats <- mclapply(seq_along(otherNames), function(x) {
  as.matrix(read.table(paste0(otherDirs[x],
                              otherNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(otherNames))

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
log2ChIPmats <- mclapply(seq_along(ChIPmats), function(x) {
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
    log2((ChIPmats[[x]]+1)/(controlmats[[1]]+1))
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIPmats[[x]]+1)/(controlmats[[2]]+1))
  }
}, mc.cores = length(ChIPmats))

for(x in seq_along(log2ChIPmats)) {
  attr(log2ChIPmats[[x]], "upstream_index") = 1:(upstream/binSize)
  attr(log2ChIPmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(log2ChIPmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(log2ChIPmats[[x]], "extend") = c(upstream, downstream)
  attr(log2ChIPmats[[x]], "smooth") = FALSE
  attr(log2ChIPmats[[x]], "signal_name") = ChIPNamesPlot[x]
  attr(log2ChIPmats[[x]], "target_name") = featureName
  attr(log2ChIPmats[[x]], "target_is_single_point") = FALSE
  attr(log2ChIPmats[[x]], "background") = 0
  attr(log2ChIPmats[[x]], "signal_is_categorical") = FALSE
  class(log2ChIPmats[[x]]) = c("normalizedMatrix", "matrix")
}

for(x in seq_along(othermats)) {
  attr(othermats[[x]], "upstream_index") = 1:(upstream/binSize)
  attr(othermats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(othermats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(othermats[[x]], "extend") = c(upstream, downstream)
  attr(othermats[[x]], "smooth") = FALSE
  attr(othermats[[x]], "signal_name") = otherNamesPlot[x]
  attr(othermats[[x]], "target_name") = featureName
  attr(othermats[[x]], "target_is_single_point") = FALSE
  attr(othermats[[x]], "background") = 0
  attr(othermats[[x]], "signal_is_categorical") = FALSE
  class(othermats[[x]]) = c("normalizedMatrix", "matrix")
}

for(x in seq_along(controlmats)) {
  attr(controlmats[[x]], "upstream_index") = 1:(upstream/binSize)
  attr(controlmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(controlmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(controlmats[[x]], "extend") = c(upstream, downstream)
  attr(controlmats[[x]], "smooth") = FALSE
  attr(controlmats[[x]], "signal_name") = controlNamesPlot[x]
  attr(controlmats[[x]], "target_name") = featureName
  attr(controlmats[[x]], "target_is_single_point") = FALSE
  attr(controlmats[[x]], "background") = 0
  attr(controlmats[[x]], "signal_is_categorical") = FALSE
  class(controlmats[[x]]) = c("normalizedMatrix", "matrix")
}

if(grepl("genes", featureName)) {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Heatmap plotting function
# Note that for plotting heatmaps for individual datasets in separate PDFs,
# must edit this function - print(EnrichedHeatmap(...))
featureHeatmap <- function(mat,
                           col_fun,
                           colour,
                           datName) {
  EnrichedHeatmap(mat = mat,
                  col = col_fun,
                  column_title = datName,
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = colour,
                                                                                        lwd = 2),
                                                                              yaxis_side = "left",
                                                                              yaxis_facing = "left",
                                                                              yaxis_gp = gpar(fontsize = 10),
                                                                              pos_line_gp = gpar(col = "black",
                                                                                                 lty = 2,
                                                                                                 lwd = 2))),
                  top_annotation_height = unit(2, "cm"),
                  width = unit(6, "cm"),
                  name = datName,
                  heatmap_legend_param = list(title = datName,
                                              title_position = "topcenter",
                                              title_gp = gpar(font = 2, fontsize = 12),
                                              legend_direction = "horizontal",
                                              labels_gp = gpar(fontsize = 10)),
                  axis_name = c(paste0("-", flankNamePlot),
                                featureStartLab, featureEndLab,
                                paste0("+", flankNamePlot)),
                  axis_name_gp = gpar(fontsize = 12),
                  border = FALSE,
                  pos_line_gp = gpar(col = "white", lty = 2, lwd = 2),
                  # If converting into png with pdfTotiffTopng.sh,
                  # set use_raster to FALSE
                  #use_raster = FALSE)
                  use_raster = TRUE, raster_device = "png", raster_quality = 10)
}

# Define heatmap colours
rich8to6equal <- c("#0000CB", "#0081FF", "#87CEFA", "#FDEE02", "#FFAB00", "#FF3300")
quantileColours <- c("darkorange1", "green2", "purple3", "deepskyblue")

# Create quantile colour block "heatmap"
quantileBlockhtmp <-   Heatmap(km$quantile,
                              col = structure(quantileColours,
                                              names = paste0("Cluster ", 1:quantiles)),
                              show_row_names = FALSE, show_heatmap_legend = FALSE,
                              width = unit(3, "mm"), name = "") 
# Plot together
log2ChIPhtmpList <- mclapply(seq_along(ChIPNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(log2ChIPmats[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(mat = log2ChIPmats[[x]],
                 col_fun = ChIP_col_fun,
                 colour = quantileColours,
                 datName = ChIPNamesPlot[x])
}, mc.cores = length(log2ChIPmats))
otherhtmpList <- mclapply(seq_along(otherNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(othermats[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(mat = othermats[[x]],
                 col_fun = ChIP_col_fun,
                 colour = quantileColours,
                 datName = otherNamesPlot[x])
}, mc.cores = length(othermats))
controlhtmpList <- mclapply(seq_along(controlNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(controlmats[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(mat = controlmats[[x]],
                 col_fun = ChIP_col_fun,
                 colour = quantileColours,
                 datName = controlNamesPlot[x])
}, mc.cores = length(controlmats))

htmpList <- c(quantileBlockhtmp, 
              log2ChIPhtmpList,
              otherhtmpList,
              controlhtmpList[[1]])

htmps <- NULL
for(x in 1:length(htmpList)) {
  htmps <- htmps + htmpList[[x]]
}
pdf(paste0(plotDir, "log2ChIPcontrol_around_", featureName,
           "_heatmaps_quantileed_by_log2_", libName, "_control_in_", region, ".pdf"),
    width = 3*length(htmpList),
    height = 8)
draw(htmps,
     split = km$quantile,
     row_order = row_order,
     heatmap_legend_side = "bottom",
     gap = unit(c(1, rep(14, length(htmpList)-1)), "mm")
    )
dev.off()

## ChIP
#ChIP_col_fun <- colorRamp2(quantile(log2ChIPmat,
#                                    c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                    na.rm = T),
#                           rich8to6equal)
#log2ChIPhtmp <- featureHeatmap(mat = log2ChIPmat,
#                               col_fun = ChIP_col_fun,
#                               colour = c("darkorange1", "green2", "purple3", "deepskyblue"),
#                               datName = "ASY1",
#                               rowSplit = km$quantile)
#pdf(paste0(plotDir, "log2ChIPcontrol_around_", featureName,
#           "_heatmaps_quantileed_by_log2_", libName, "_control_in_", region, ".pdf"),
#    width = 3,
#    height = 8)
#draw(log2ChIPhtmp,
#     split = km$quantile,
#     row_order = row_order,
#     heatmap_legend_side = "bottom",
#     gap = unit(c(2, 14), "mm")
#    )
#dev.off()
