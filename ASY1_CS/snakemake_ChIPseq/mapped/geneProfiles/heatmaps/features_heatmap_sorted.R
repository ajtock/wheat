#!/applications/R/R-3.4.0/bin/Rscript

# Plot heatmaps of features sorted by coverage levels in promoters,
# within feature bodies, or in terminators

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_heatmap_sorted.R ASY1_CS_Rep1_ChIP both 107891 genes 3500 2000 2kb '2 kb' 20 20bp promoters

libName <- "ASY1_CS_Rep1_ChIP"
align <- "both"
featureNumber <- 107891
featureName <- "genes"
bodyLength <- 3500
upstream <- 2000
downstream <- 2000
flankName <- "2kb"
flankNamePlot <- "2 kb"
binSize <- 20
binName <- "20bp"
region <- "promoters"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
align <- args[2]
featureNumber <- as.numeric(args[3])
featureName <- args[4]
bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[7])
flankName <- args[8]
flankNamePlot <- args[9]
binSize <- as.numeric(args[10])
binName <- args[11]
region <- args[12]

library(EnrichedHeatmap)
library(parallel)
library(circlize)
library(RColorBrewer)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

plotDir <- "./plots/"
regionPlotDir <- paste0(plotDir, region, "/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
system(paste0("[ -d ", regionPlotDir, " ] || mkdir ", regionPlotDir))
plotDir <- regionPlotDir
 
# Load matrix and extract region for sorting of features
mat1 <- as.matrix(read.table(paste0("../matrices/",
                                    libName,
                                    "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                    featureName, "_matrix_bin", binName,
                                    "_flank", flankName, ".tab"),
                             header = F, skip = 3))
if( region == "promoters" ) {
  mat1Region <- mat1[,(((upstream-500)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  mat1Region <- mat1[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(500/binSize))]
} else if ( region == "bodies" ) {
  mat1Region <- mat1[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else {
  print("The region name provided does not match 'promoters', 'terminators', or 'bodies'")
}
mat1RegionRowMeans <- rowMeans(mat1Region, na.rm = T)
mat1RegionRowMeansSorted <- sort.int(mat1RegionRowMeans,
                                     decreasing = T,
                                     index.return = T,
                                     na.last = T)
mat1RegionRowSums <- rowSums(mat1Region, na.rm = T)
mat1RegionRowSumsSorted <- sort.int(mat1RegionRowSums,
                                    decreasing = T,
                                    index.return = T,
                                    na.last = T)
mat1RegionRowMedians <- apply(X = mat1Region,
                              MARGIN = 1,
                              FUN = median)
mat1RegionRowMediansSorted <- sort.int(mat1RegionRowMedians,
                                       decreasing = T,
                                       index.return = T,
                                       na.last = T)
mat1Sorted <- mat1[sort.int(mat1RegionRowMeans,
                            decreasing = T,
                            index.return = T,
                            na.last = T)$ix,]

# Load feature matrices for each chromatin dataset, log2-transform,
# and sort by decreasing mat1RegionRowMeans

ChIPNames1 <- c(
                "
                "H3K9me2_Rep1_ChIP",
                "H3K4me3_Rep1_ChIP",
                "CENH3_ChIP_SRR1686799"
               )
ChIPNames2 <- c(
                "H3K4me3_ChIP_SRR6350668",
                "H3K9ac_ChIP_SRR6350667",
                "H3K27me3_ChIP_SRR6350666",
                "H3K36me3_ChIP_SRR6350670"
               )
controlNames <- c(
                  "MNase_Rep1",
                  "H3_input_SRR6350669"
                 )
ChIPNamesPlot <- c(
                  "H3K9me2",
                  "H3K4me3",
                  "CENH3",
                  "H3K4me3 (IWGSC)",
                  "H3K9ac",
                  "H3K27me3",
                  "H3K36me3",
                  "MNase"
                 )
ChIPNames <- c(
               "H3K9me2_Rep1_ChIP",
               "H3K4me3_Rep1_ChIP",
               "CENH3_ChIP_SRR1686799",
               "H3K4me3_ChIP_SRR6350668",
               "H3K9ac_ChIP_SRR6350667",
               "H3K27me3_ChIP_SRR6350666",
               "H3K36me3_ChIP_SRR6350670",
               "MNase_Rep1"
              )
ChIPDirs1 <- c(
               "/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/CENH3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"
              )
ChIPDirs2 <- c(
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K9ac/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K36me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"
              )
controlDirs <- c(
                 "/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
                 "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/input/snakemake_ChIPseq/mapped/geneProfiles/matrices/"             
                )

ChIPmats1 <- mclapply(seq_along(ChIPNames1), function(x) {
  as.matrix(read.table(paste0(ChIPDirs1[x],
                              ChIPNames1[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames1)) 

ChIPmats2 <- mclapply(seq_along(ChIPNames2), function(x) {
  as.matrix(read.table(paste0(ChIPDirs2[x],
                              ChIPNames2[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames2)) 

controlmats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames)) 

log2ChIPmats1 <- mclapply(seq_along(ChIPmats1), function(x) {
  log2((ChIPmats1[[x]]+1)/(controlmats[[1]]+1))
}, mc.cores = length(ChIPmats1))
 
log2ChIPmats2 <- mclapply(seq_along(ChIPmats2), function(x) {
  log2((ChIPmats2[[x]]+1)/(controlmats[[2]]+1))
}, mc.cores = length(ChIPmats2))

log2ChIPmats <- c(log2ChIPmats1, log2ChIPmats2, list(controlmats[[1]]))

log2ChIPmatsSorted <- mclapply(seq_along(log2ChIPmats), function(x) {
  log2ChIPmats[[x]][sort.int(mat1RegionRowMeans,
                             decreasing = T,
                             index.return = T,
                             na.last = T)$ix,]
}, mc.cores = length(log2ChIPmats))

for(x in seq_along(log2ChIPmatsSorted)) {
  attr(log2ChIPmatsSorted[[x]], "upstream_index") = 1:(upstream/binSize)
  attr(log2ChIPmatsSorted[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(log2ChIPmatsSorted[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(log2ChIPmatsSorted[[x]], "extend") = c(upstream, downstream)
  attr(log2ChIPmatsSorted[[x]], "smooth") = FALSE
  attr(log2ChIPmatsSorted[[x]], "signal_name") = ChIPNamesPlot[x]
  attr(log2ChIPmatsSorted[[x]], "target_name") = featureName
  attr(log2ChIPmatsSorted[[x]], "target_is_single_point") = FALSE
  attr(log2ChIPmatsSorted[[x]], "background") = 0
  attr(log2ChIPmatsSorted[[x]], "signal_is_categorical") = FALSE
  class(log2ChIPmatsSorted[[x]]) = c("normalizedMatrix", "matrix")
}

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Heatmap plotting function
geneHeatmap <- function(matSorted,
                        col_fun,
                        datName,
                        rowOrder
                       ) {
  ht <- print(EnrichedHeatmap(mat = matSorted,
                              top_annotation = NULL,
                              col = col_fun,
                              name = "Coverage",
                              row_order = rowOrder,
                              column_title = datName,
                              axis_name = c(paste0("-", flankNamePlot),
                                            featureStartLab, featureEndLab,
                                            paste0("+", flankNamePlot)),
                              border = FALSE,
                              pos_line_gp = gpar(col = "white", lty = 2, lwd = 2),
                              use_raster = TRUE))
  
}

# Plot sRNA data
rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")
foreach(x = 1:5) %dopar% {
  sRNA_col_fun <- colorRamp2(quantile(c(sRNAmatsSorted[[x]]),
                                      c(0.998, 0.9982, 0.9984, 0.9986, 0.9988, 0.999),
                                      na.rm = T),
                             rich8to6equal)
  pdf(paste0(plotDir, sRNAsizes[x], "nt_sRNAs_around_top", featureNumber, featureName,
             "_heatmap_ordered_by_",
             libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  geneHeatmap(matSorted = sRNAmatsSorted[[x]][1:featureNumber,],
              col_fun = sRNA_col_fun,
              datName = paste0(sRNAsizes[x], "-nt sRNAs"),
              rowOrder = c(1:dim(sRNAmatsSorted[[x]][1:featureNumber,])[1])
             )
  dev.off()
}

# 33-nt sRNAs
  sRNA_col_fun <- colorRamp2(quantile(c(sRNAmatsSorted[[6]]),
                                      c(0.999, 0.99902, 0.99904, 0.99906, 0.99908, 0.9991),
                                      na.rm = T),
                             rich8to6equal)
  pdf(paste0(plotDir, sRNAsizes[6], "nt_sRNAs_around_top", featureNumber, featureName,
             "_heatmap_ordered_by_",
             libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  geneHeatmap(matSorted = sRNAmatsSorted[[6]][1:featureNumber,],
              col_fun = sRNA_col_fun,
              datName = paste0(sRNAsizes[6], "-nt sRNAs"),
              rowOrder = c(1:dim(sRNAmatsSorted[[6]][1:featureNumber,])[1])
             )
  dev.off()

# 34-nt sRNAs
  sRNA_col_fun <- colorRamp2(quantile(c(sRNAmatsSorted[[7]]),
                                      c(0.999, 0.9992, 0.9994, 0.9996, 0.99985, 0.9999),
                                      na.rm = T),
                             rich8to6equal)
  pdf(paste0(plotDir, sRNAsizes[7], "nt_sRNAs_around_top", featureNumber, featureName,
             "_heatmap_ordered_by_",
             libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  geneHeatmap(matSorted = sRNAmatsSorted[[7]][1:featureNumber,],
              col_fun = sRNA_col_fun,
              datName = paste0(sRNAsizes[7], "-nt sRNAs"),
              rowOrder = c(1:dim(sRNAmatsSorted[[7]][1:featureNumber,])[1])
             )
  dev.off()

# Plot ChIP data
foreach(x = 1:7) %dopar% {
  ChIP_col_fun <- colorRamp2(quantile(c(log2ChIPmatsSorted[[x]]),
                                      c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  pdf(paste0(plotDir, ChIPNames[x], "_around_top", featureNumber, featureName,
             "_heatmap_ordered_by_",
             libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  geneHeatmap(matSorted = log2ChIPmatsSorted[[x]][1:featureNumber,],
              col_fun = ChIP_col_fun,
              datName = ChIPNamesPlot[x],
              rowOrder = c(1:dim(log2ChIPmatsSorted[[x]][1:featureNumber,])[1])
             )
  dev.off()
}

# MNase
  ChIP_col_fun <- colorRamp2(quantile(c(log2ChIPmatsSorted[[8]]),
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  pdf(paste0(plotDir, ChIPNames[8], "_around_top", featureNumber, featureName,
             "_heatmap_ordered_by_",
             libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  geneHeatmap(matSorted = log2ChIPmatsSorted[[8]][1:featureNumber,],
              col_fun = ChIP_col_fun,
              datName = ChIPNamesPlot[8],
              rowOrder = c(1:dim(log2ChIPmatsSorted[[8]][1:featureNumber,])[1])
             )
  dev.off()
