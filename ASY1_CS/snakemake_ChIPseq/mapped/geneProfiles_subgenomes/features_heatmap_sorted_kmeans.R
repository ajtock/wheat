#!/applications/R/R-3.4.0/bin/Rscript

# Cluster genes by log2(ChIP/control)

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_heatmap_sorted_kmeans.R ASY1_CS_Rep1_ChIP ASY1_CS both genes_in_Agenome_genomewide 3500 2000 2kb '2 kb' 20 20bp promoters

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#align <- "both"
#featureName <- "genes_in_Agenome_genomewide"
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
featureName <- args[4]
bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
flankNamePlot <- args[8]
binSize <- as.numeric(args[9])
binName <- args[10]
region <- args[11]

library(EnrichedHeatmap)
library(circlize)
library(RColorBrewer)
library(gridExtra)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

outDir <- paste0("clusters_by_log2_", libName,
                 "_control_in_", region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load ChIP matrix
mat1 <- as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/",
                                    dirName,
                                    "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                    libName,
                                    "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                    featureName, "_matrix_bin", binName,
                                    "_flank", flankName, ".tab"),
                             header = F, skip = 3))

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
controlmats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

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

# Extract region for clustering of features (adjust promoter/terminator size as necessary)
if( region == "promoters" ) {
  log2ChIPmatRegion <- log2ChIPmat[,(((upstream-500)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  log2ChIPmatRegion <- log2ChIPmat[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(500/binSize))]
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
# Replace NAs in log2ChIPmatRegion with 0
log2ChIPmatRegion[which(is.na(log2ChIPmatRegion))] <- 0

# Determine number of clusters to be used for k-means clustering
# by generating a scree ("elbow") plot of the ratio of the
# within-cluster sum of squares (WSS) to the total sum of squares (TSS)
# for k clusters

# First run garbage collection to free up memory
gc()
# Set seed for reproducible clusters
set.seed(938402845)

# Initialise ratio_ss as a vector of 0s to be replaced with the ratio
# of the within-cluster sum of squares to the total sum of squares
# for each number of clusters
kMax <- 10
ratio_ss <- rep(0, kMax)

# Apply k-means clustering to log2ChIPmat for k in 1:15
# and obtain the ratio of the within-cluster sum of squares (WSS) to
# the total sum of squares (TSS) for each number of clusters
for(k in 1:kMax) {
  print(k)
  km <- kmeans(x = log2ChIPmatRegion,
               centers = k,
               iter.max = 10,
               nstart = 10)
  ratio_ss[k] <- km$tot.withinss / km$totss
}
#19: Quick-TRANSfer stage steps exceeded maximum (= 1767250)

# Make a scree ("elbow") plot
pdf(paste0(outDir, "screePlot_k_clusters_by_log2_",
           libName, "_control_in_",
           region, "_of_", featureName, ".pdf"))
plot(x = 1:kMax, y = ratio_ss,
     type = "b",
     xlab = "k", ylab = "WSS / TSS")
abline(h = 0.2, lty = 2, col = "red")
dev.off()

# First run garbage collection to free up memory
gc()
# Set seed for reproducible clusters
set.seed(938402845)

# 4 clusters
kDef <- 4
km <- kmeans(x = log2ChIPmatRegion,
             centers = kDef,
             iter.max = 10,
             nstart = 10)
# Adjust the cluster numbers so that cluster 1 has highest levels
x <- tapply(X = rowMeans(log2ChIPmatRegion),
            INDEX = km$cluster,
            FUN = mean)
od <- order(structure(order(x, decreasing = TRUE),
                      names = names(x)),
            decreasing = F)
km$cluster <- od[km$cluster]
km$cluster <- paste0("Cluster ", km$cluster)

## Calculate Dunn's index to assess compactness and separation of clusters
## "The Dunn Index is the ratio of the smallest distance between observations
## not in the same cluster to the largest intra-cluster distance.
## The Dunn Index has a value between zero and infinity, and should be maximized."
## (NOTE: takes a long time)
#library(clValid)
#km_dunn <- dunn(clusters = km$cluster, Data = log2ChIPmatRegion)

# Order genes in each cluster by decreasing log2ChIPmatRegion levels
# to define "row_order" for heatmaps
combineRowOrders <- function(cluster_bool_list) {
  do.call("c", lapply(cluster_bool_list, function(x) {
    cluster_log2ChIPmatRegionRowMeans <- rowMeans(log2ChIPmatRegion[x,], na.rm = T)
    which(x)[order(cluster_log2ChIPmatRegionRowMeans, decreasing = T)]
  }))
}
row_order <- combineRowOrders(cluster_bool_list =
  lapply(seq_along(1:kDef), function(k) { 
    km$cluster == paste0("Cluster ", k)
  })
)
# Order gene IDs in each cluster by decreasing log2ChIPmatRegion levels
# for use in GO term enrichment analysis
listCombineRowOrders <- function(cluster_bool_list) {
  do.call(list, lapply(cluster_bool_list, function(x) {
    cluster_log2ChIPmatRegionRowMeans <- rowMeans(log2ChIPmatRegion[x,], na.rm = T)
    which(x)[order(cluster_log2ChIPmatRegionRowMeans, decreasing = T)]
  }))
}
featureIndicesList <- listCombineRowOrders(cluster_bool_list =
  lapply(seq_along(1:kDef), function(k) {
    km$cluster == paste0("Cluster ", k)
  })
)
# Alternatively, with original ordering:
## Get feature indices for each cluster
#featureIndicesList <- lapply(seq_along(1:kDef), function(k) {
#  which(km$cluster == paste0("Cluster ", k))
#})

# Load features 
features <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                              substring(featureName, first = 10), ".gff3"),
                       header = F)
# Separate into clusters
featuresClusterList <- lapply(seq_along(1:kDef), function(k) {
  features[featureIndicesList[[k]],]
})
featureIDsClusterList <- lapply(seq_along(1:kDef), function(k) {
  sub(pattern = "\\.\\d+", replacement = "",
      x = as.vector(features[featureIndicesList[[k]],]$V9))
})
sapply(seq_along(featureIDsClusterList), function(k) {
  write.table(featureIDsClusterList[[k]],
              file = paste0(outDir, "cluster", as.character(k), "_of_", as.character(kDef),
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
clusterColours <- c("darkorange1", "green2", "purple3", "deepskyblue")

# Create cluster colour block "heatmap"
clusterBlockhtmp <-   Heatmap(km$cluster,
                              col = structure(clusterColours,
                                              names = paste0("Cluster ", 1:kDef)),
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
                 colour = clusterColours,
                 datName = ChIPNamesPlot[x])
}, mc.cores = length(log2ChIPmats))
otherhtmpList <- mclapply(seq_along(otherNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(othermats[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(mat = othermats[[x]],
                 col_fun = ChIP_col_fun,
                 colour = clusterColours,
                 datName = otherNamesPlot[x])
}, mc.cores = length(othermats))
controlhtmpList <- mclapply(seq_along(controlNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(controlmats[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(mat = controlmats[[x]],
                 col_fun = ChIP_col_fun,
                 colour = clusterColours,
                 datName = controlNamesPlot[x])
}, mc.cores = length(controlmats))

htmpList <- c(clusterBlockhtmp, 
              log2ChIPhtmpList,
              otherhtmpList,
              controlhtmpList[[1]])

htmps <- NULL
for(x in 1:length(htmpList)) {
  htmps <- htmps + htmpList[[x]]
}
pdf(paste0(plotDir, "log2ChIPcontrol_around_", featureName,
           "_heatmaps_clustered_by_log2_", libName, "_control_in_", region, ".pdf"),
    width = 3*length(htmpList),
    height = 8)
draw(htmps,
     split = km$cluster,
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
#                               rowSplit = km$cluster)
#pdf(paste0(plotDir, "log2ChIPcontrol_around_", featureName,
#           "_heatmaps_clustered_by_log2_", libName, "_control_in_", region, ".pdf"),
#    width = 3,
#    height = 8)
#draw(log2ChIPhtmp,
#     split = km$cluster,
#     row_order = row_order,
#     heatmap_legend_side = "bottom",
#     gap = unit(c(2, 14), "mm")
#    )
#dev.off()
