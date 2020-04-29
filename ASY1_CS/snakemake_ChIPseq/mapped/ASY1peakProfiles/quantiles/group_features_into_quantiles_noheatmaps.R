#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 28.04.2020

#
# Divide features into quantiles based on mean recombination rate (cM/Mb)
# Extract and save feature IDs for each quantile for further analyses
# (e.g., average + 95% CI profile plotting).
# Calculate mean winName-scaled recombination rate (cM/Mb) from
# the promoter to the terminator of each feature.
#

# Usage:
# /applications/R/R-3.4.0/bin/Rscript group_features_into_quantiles_noheatmaps.R ASY1_CS_Rep1_ChIP ASY1_CS 'Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin' 4

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#featureName <- unlist(strsplit("Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin",
#                               split = ","))
#quantiles <- 4

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
quantiles <- as.numeric(args[4])

library(GenomicRanges)
library(dplyr)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Load features in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
features <- lapply(seq_along(featureName), function(y) {
  tmp <- read.table(paste0("/home/ajt200/analysis/wheat/", dirName,
                           "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           libName,
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           featureName[y], ".bed"),
                    header = F)
  data.frame(tmp,
             V7 = paste0(featureName[y], "_", tmp$V4),
             stringsAsFactors = F) 
})
if(length(featureName) > 1) {
  features <- do.call(rbind, features)
} else {
  features <- features[[1]]
}
colnames(features) <- c("chr", "start", "end", "name", "score", "strand", "featureID")
featuresGR <- GRanges(seqnames = features$chr,
                      ranges = IRanges(start = features$start+1,
                                       end = features$end),
                      strand = features$strand,
                      featureID = features$featureID)
# Extend feature boundaries to include 1000 bp upstream
# and downstream for calculation of mean cM/Mb around loci
featuresGR_ext <- GRanges(seqnames = seqnames(featuresGR),
                          ranges = IRanges(start = start(featuresGR)-1000,
                                           end = end(featuresGR)+1000),
                          strand = strand(featuresGR),
                          featureID = featuresGR$featureID)
print(featuresGR_ext)

# Load ranLocs in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLocs <- lapply(seq_along(featureName), function(y) {
  tmp <- read.table(paste0("/home/ajt200/analysis/wheat/", dirName,
                           "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           libName,
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           featureName[y], "_randomLoci.bed"),
                    header = F)
  data.frame(tmp,
             V7 = paste0(featureName[y], "_", tmp$V4),
             stringsAsFactors = F) 
})
if(length(featureName) > 1) {
  ranLocs <- do.call(rbind, ranLocs)
} else {
  ranLocs <- ranLocs[[1]]
}
colnames(ranLocs) <- c("chr", "start", "end", "name", "score", "strand", "ranLocID")
ranLocsGR <- GRanges(seqnames = ranLocs$chr,
                     ranges = IRanges(start = ranLocs$start+1,
                                      end = ranLocs$end),
                     strand = ranLocs$strand,
                     ranLocID = ranLocs$ranLocID)
# Extend ranLoc boundaries to include 1000 bp upstream
# and downstream for calculation of mean cM/Mb around loci
ranLocsGR_ext <- GRanges(seqnames = seqnames(ranLocsGR),
                         ranges = IRanges(start = start(ranLocsGR)-1000,
                                          end = end(ranLocsGR)+1000),
                         strand = strand(ranLocsGR),
                         ranLocID = ranLocsGR$ranLocID)
print(ranLocsGR_ext)

# Convert windowed recombination rate into GRanges
cMMb <- read.table(paste0(
                   "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_recombination_rate_analysis/",
                   "iwgsc_refseqv1.0_recombination_rate.txt"),
                   header = T)
cMMbGR <- GRanges(seqnames = cMMb$chromosome,
                  ranges = IRanges(start = cMMb$intervalStart,
                                   end = cMMb$intervalEnd),
                  strand = "*",
                  cMMb = cMMb$recombinationRate)

# Obtain cMMb values for each feature extended by 1000 bp on each side
# Where features overlap more than one winName window, calculate mean cMMb
feature_cMMb_overlaps <- findOverlaps(query = featuresGR_ext,
                                      subject = cMMbGR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
feature_cMMb_overlapsList <- lapply(seq_along(featuresGR_ext), function(x) {
  subjectHits(feature_cMMb_overlaps)[queryHits(feature_cMMb_overlaps) == x]
})
feature_cMMb <- sapply(feature_cMMb_overlapsList,
                       function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
featuresGR <- GRanges(featuresGR,
                      featureID = featuresGR$featureID,
                      cMMb = feature_cMMb)

# Obtain cMMb values for each ranLoc extended by 1000 bp on each side
# Where ranLocs overlap more than one winName window, calculate mean cMMb
ranLoc_cMMb_overlaps <- findOverlaps(query = ranLocsGR_ext,
                                     subject = cMMbGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
ranLoc_cMMb_overlapsList <- lapply(seq_along(ranLocsGR_ext), function(x) {
  subjectHits(ranLoc_cMMb_overlaps)[queryHits(ranLoc_cMMb_overlaps) == x]
})
ranLoc_cMMb <- sapply(ranLoc_cMMb_overlapsList,
                      function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
ranLocsGR <- GRanges(ranLocsGR,
                     ranLocID = ranLocsGR$ranLocID,
                     cMMb = ranLoc_cMMb)


# Define orderingFactor to be used for grouping features into quantiles
orderingFactor <- "cMMb"
outDir <- paste0("quantiles_by_", orderingFactor, "/")
plotDir_list <- lapply(seq_along(outDir), function(w) {
  paste0(outDir[w], "plots/")
})

sapply(seq_along(outDir), function(w) {
 system(paste0("[ -d ", outDir[w], " ] || mkdir ", outDir[w]))
})
sapply(seq_along(outDir), function(w) {
  system(paste0("[ -d ", plotDir_list[[w]], " ] || mkdir ", plotDir_list[[w]]))
})

# For each population, divide features into quantiles based on decreasing orderingFactor
# all subgenomes
featuresDF <- data.frame(featuresGR,
                         quantile = as.character(""),
                         stringsAsFactors = F)
mclapply(seq_along(orderingFactor), function(w) {
  print(orderingFactor[w])
  # Assign 0s to NA values only for coverage data
  if(grepl("_in_", orderingFactor[w])) {
    featuresDF[,which(colnames(featuresDF) == orderingFactor[w])][which(is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]))] <- 0
  }
  if(grepl("HudsonRM", orderingFactor[w])) {
    quantiles <- 2
  } else {
    quantiles <- 4
  }
  quantilesStats <- data.frame()
  for(k in 1:quantiles) {
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
    if(k < quantiles) {
      featuresDF[ !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >  1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    } else {
    # Final quantile should span 0 to, e.g., 0.25 proportions of features
      featuresDF[ !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                  percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >= 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
    }
    write.table(featuresDF[featuresDF$quantile == paste0("Quantile ", k),],
                file = paste0(outDir_list[[w]][x],
                              "quantile", k, "_of_", quantiles,
                              "_by_", orderingFactor[w],
                              "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    stats <- data.frame(quantile = as.integer(k),
                        n = as.integer(dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]),
                        mean_width = as.integer(round(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                        total_width = as.integer(sum(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                        mean_orderingFactor = as.numeric(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),][,which(colnames(featuresDF) == orderingFactor[w])], na.rm = T)))
    quantilesStats <- rbind(quantilesStats, stats)
  }
  write.table(quantilesStats,
              file = paste0(outDir_list[[w]][x],
                            "summary_", quantiles, "quantiles_by_", orderingFactor[w], "_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(featuresDF,
              file = paste0(outDir_list[[w]][x],
                            "features_", quantiles, "quantiles",
                            "_by_", orderingFactor[w], "_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # Divide ranLocs into quantiles based on feature quantile indices
  ranLocsDF <- data.frame(ranLocsGR,
                          random = as.character(""),
                          stringsAsFactors = F)
  # Get row indices for each feature quantile
  quantileIndices <- lapply(1:quantiles, function(k) {
    which(featuresDF$quantile == paste0("Quantile ", k))
  })
  for(k in 1:quantiles) {
    ranLocsDF[quantileIndices[[k]],]$random <- paste0("Random ", k)
  }
  write.table(ranLocsDF,
              file = paste0(outDir_list[[w]][x],
                            "features_", quantiles, "quantiles",
                            "_by_", orderingFactor[w], "_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), "_", pop_name[x], "_ranLocs.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}, mc.cores = length(orderingFactor), mc.preschedule = F)



# Divide features into quantiles based on decreasing log2ChIPmatRegionRowMeans
featuresDF <- data.frame(featuresGR,
                         quantile = as.character(""),
                         stringsAsFactors = F)
featuresDF$log2ChIPmatRegionRowMeans[which(is.na(featuresDF$log2ChIPmatRegionRowMeans))] <- 0
quantilesStats <- data.frame()
for(k in 1:quantiles) {
  if(k < quantiles) {
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
    featuresDF[ percent_rank(featuresDF$log2ChIPmatRegionRowMeans) <= 1-((k-1)/quantiles) &
                percent_rank(featuresDF$log2ChIPmatRegionRowMeans) >  1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of features
    featuresDF[ percent_rank(featuresDF$log2ChIPmatRegionRowMeans) <= 1-((k-1)/quantiles) &
                percent_rank(featuresDF$log2ChIPmatRegionRowMeans) >= 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  }
  write.table(featuresDF[featuresDF$quantile == paste0("Quantile ", k),],
              file = paste0(outDir,
                            "quantile", k, "_of_", quantiles,
                            "_by_log2_", libName, "_control_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  stats <- data.frame(quantile = as.integer(k),
                      n = as.integer(dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]),
                      mean_width = as.integer(round(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                      total_width = as.integer(sum(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                      mean_log2ChIPmatRegionRowMeans = as.numeric(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$log2ChIPmatRegionRowMeans, na.rm = T)))
  quantilesStats <- rbind(quantilesStats, stats)
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
write.table(featuresDF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control_in_",
                          region, "_of_",
                          substring(featureName[1][1], first = 1, last = 5), "_in_",
                          paste0(substring(featureName, first = 10, last = 16),
                                 collapse = "_"), "_",
                          substring(featureName[1][1], first = 18), ".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Divide ranLocs into quantiles based on feature quantile indices
ranLocsDF <- data.frame(ranLocsGR,
                        random = as.character(""),
                        stringsAsFactors = F)
# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$quantile == paste0("Quantile ", k))
})
for(k in 1:quantiles) {
  ranLocsDF[quantileIndices[[k]],]$random <- paste0("Random ", k)
}
write.table(ranLocsDF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control_in_",
                          region, "_of_",
                          substring(featureName[1][1], first = 1, last = 5), "_in_",
                          paste0(substring(featureName, first = 10, last = 16),
                                 collapse = "_"), "_",
                          substring(featureName[1][1], first = 18), "_ranLocs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Order features in each quantile by decreasing log2ChIPmatRegion levels
# to define "row_order" for heatmaps
combineRowOrders <- function(quantile_bool_list) {
  do.call("c", lapply(quantile_bool_list, function(x) {
    quantile_log2ChIPmatRegionRowMeans <- rowMeans(log2ChIPmatRegion[x,], na.rm = T)
    quantile_log2ChIPmatRegionRowMeans[which(is.na(quantile_log2ChIPmatRegionRowMeans))] <- 0
    which(x)[order(quantile_log2ChIPmatRegionRowMeans, decreasing = T)]
  }))
}
row_order <- combineRowOrders(quantile_bool_list =
  lapply(seq_along(1:quantiles), function(k) { 
    featuresDF$quantile == paste0("Quantile ", k)
  })
)
# Confirm row_order is as would be obtained by alternative method
# Note that this alternative 
stopifnot(identical(row_order,
                    order(featuresDF$log2ChIPmatRegionRowMeans,
                          decreasing=T)))

# Order feature IDs in each quantile by decreasing log2ChIPmatRegion levels
# for use in GO term enrichment analysis
listCombineRowOrders <- function(quantile_bool_list) {
  do.call(list, lapply(quantile_bool_list, function(x) {
    quantile_log2ChIPmatRegionRowMeans <- rowMeans(log2ChIPmatRegion[x,], na.rm = T)
    quantile_log2ChIPmatRegionRowMeans[which(is.na(quantile_log2ChIPmatRegionRowMeans))] <- 0
    which(x)[order(quantile_log2ChIPmatRegionRowMeans, decreasing = T)]
  }))
}
featureIndicesList <- listCombineRowOrders(quantile_bool_list =
  lapply(seq_along(1:quantiles), function(k) {
    featuresDF$quantile == paste0("Quantile ", k)
  })
)
stopifnot(identical(row_order,
                    do.call(c, lapply(featureIndicesList,
                                      function(x) x))))
# Alternatively, with original ordering:
## Get feature indices for each quantile
#featureIndicesList <- lapply(seq_along(1:quantiles), function(k) {
#  which(featuresDF$quantile == paste0("Quantile ", k))
#})

featureIDsQuantileList <- lapply(seq_along(1:quantiles), function(k) {
  sub(pattern = "\\.\\d+", replacement = "",
      x = as.vector(featuresDF[featureIndicesList[[k]],]$featureID))
})
sapply(seq_along(featureIDsQuantileList), function(k) {
  write.table(featureIDsQuantileList[[k]],
              file = paste0(outDir,
                            "featureIDs_quantile", k, "_of_", quantiles,
                            "_by_log2_", libName, "_control_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), ".txt"),
              quote = F, row.names = F, col.names = F)
})
