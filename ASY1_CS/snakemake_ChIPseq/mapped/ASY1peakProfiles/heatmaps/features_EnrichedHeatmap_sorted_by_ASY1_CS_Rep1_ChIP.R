#!/applications/R/R-3.4.0/bin/Rscript

# Plot heatmaps of features sorted by coverage levels in promoters,
# within feature bodies, or in terminators

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_EnrichedHeatmap_sorted_by_ASY1_CS_Rep1_ChIP.R ASY1_CS_Rep1_ChIP ASY1_CS_peaks_in_Agenome_euchromatin 2000 2kb '2 kb' 20 20bp bodies

#libName <- "ASY1_CS_Rep1_ChIP"
#featureName <- "ASY1_CS_peaks_in_Agenome_euchromatin"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "bodies"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
featureName <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
region <- args[8]

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

plotDir <- "plots/"
regionPlotDir <- paste0(plotDir, region, "_by_", libName, "/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
system(paste0("[ -d ", regionPlotDir, " ] || mkdir ", regionPlotDir))
plotDir <- regionPlotDir
 
# Load feature covergae matrices for each chromatin dataset,
# log2-transform, and sort by decreasing log2mat1RegionRowMeans
libNameChIP <- c(
                 "H2AZ_Rep1_ChIP",
                 "H3K4me3_Rep1_ChIP",
                 "H3K9me2_Rep1_ChIP",
                 "H3K27me1_Rep1_ChIP",
                 "ASY1_CS_Rep1_ChIP",
                 "MNase_Rep1",
                 "H3K4me3_ChIP_SRR6350668",
                 "H3K27me3_ChIP_SRR6350666",
                 "H3K36me3_ChIP_SRR6350670",
                 "H3K9ac_ChIP_SRR6350667"
                )
libNameControl <- c(
                    "MNase_Rep1",
                    "H3_input_SRR6350669"
                   )
markChIP <- c(
              "H2AZ",
              "H3K4me3",
              "H3K9me2",
              "H3K27me1",
              "ASY1_CS",
              "MNase",
              "H3K4me3",
              "H3K27me3",
              "H3K36me3",
              "H3K9ac"
             )
markControl <- c(
                 "MNase",
                 "input"
                )
libNamePlot <- c(
                 "H2A.Z",
                 "H3K4me3",
                 "H3K9me2",
                 "H3K27me1",
                 "ASY1 (CS)",
                 "MNase",
                 "H3K4me3",
                 "H3K27me3",
                 "H3K36me3",
                 "H3K9ac"
                )
## ChIP
covDirChIP <- sapply(seq_along(libNameChIP), function(x) {
  if(libNameChIP[x] %in% c("H3K4me3_ChIP_SRR6350668",
                           "H3K27me3_ChIP_SRR6350666",
                           "H3K36me3_ChIP_SRR6350670",
                           "H3K9ac_ChIP_SRR6350667",
                           "CENH3_ChIP_SRR1686799")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           markChIP[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else { 
    paste0("/home/ajt200/analysis/wheat/",
           markChIP[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  }
})

covMatChIP <- mclapply(seq_along(libNameChIP), function(x) {
  as.matrix(read.table(paste0(covDirChIP[x],
                              libNameChIP[x],
                              "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(libNameChIP))

## Control
covDirControl <- sapply(seq_along(libNameControl), function(x) {
  if(libNameControl[x] == "H3_input_SRR6350669") {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           "input/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else if(libNameControl[x] == "MNase_Rep1") {
    paste0("/home/ajt200/analysis/wheat/",
           "MNase/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else {
    if(!(libNameControl %in% c("H3_input_SRR6350669", "MNase_Rep1"))) {
      stop("libNameControl[x] is neither H3_input_SRR6350669 nor MNase_Rep1")
    }
  }
})

covMatControl <- mclapply(seq_along(libNameControl), function(x) {
  as.matrix(read.table(paste0(covDirControl[x],
                              libNameControl[x],
                              "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(libNameControl))

## Load matrix, log2-transform, and extract region for sorting of features
mat1 <- as.matrix(read.table(paste0("../matrices/",
                                    libName,
                                    "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                                    featureName, "_matrix_bin", binName,
                                    "_flank", flankName, ".tab"),
                             header = F, skip = 3))

# Calculate log2(ChIP/control) for mat1, to be used for sorting features
# by decreasing levels
if(libName %in% c("ASY1_CS_Rep1_ChIP",
                  "MNase_Rep1",
                  "H3K4me3_ChIP_SRR6350668",
                  "H3K27me3_ChIP_SRR6350666",
                  "H3K36me3_ChIP_SRR6350670",
                  "H3K9ac_ChIP_SRR6350667",
                  "CENH3_ChIP_SRR1686799")) {
  # log2(ChIP/input)
  log2mat1 <- log2((mat1+1)/(covMatControl[[2]]+1))
} else {
  # log2(ChIP/MNase)
  log2mat1 <- log2((mat1+1)/(covMatControl[[1]]+1))
}

bodyLength <- (dim(log2mat1)[2]-((upstream+downstream)/binSize))*binSize
if( region == "promoters" ) {
  log2mat1Region <- log2mat1[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  log2mat1Region <- log2mat1[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
} else if ( region == "bodies" ) {
  log2mat1Region <- log2mat1[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else {
  print("The region name provided does not match 'promoters',  'terminators', or 'bodies'")
}
log2mat1RegionRowMeans <- rowMeans(log2mat1Region, na.rm = T)
log2mat1RegionRowMeansSorted <- sort.int(log2mat1RegionRowMeans,
                                         decreasing = T,
                                         index.return = T,
                                         na.last = T)

## log2(ChIP/control) matrices
log2covMat <- mclapply(seq_along(libNameChIP), function(x) {
  if(libNameChIP[x] %in% c("ASY1_CS_Rep1_ChIP",
                           "MNase_Rep1",
                           "H3K4me3_ChIP_SRR6350668",
                           "H3K27me3_ChIP_SRR6350666",
                           "H3K36me3_ChIP_SRR6350670",
                           "H3K9ac_ChIP_SRR6350667",
                           "CENH3_ChIP_SRR1686799")) {
    # log2(ChIP/input)
    log2((covMatChIP[[x]]+1)/(covMatControl[[2]]+1))
  } else {
    # log2(ChIP/MNase)
    log2((covMatChIP[[x]]+1)/(covMatControl[[1]]+1))
  }
}, mc.cores = length(libNameChIP))

log2covMatSorted <- mclapply(seq_along(log2covMat), function(x) {
  log2covMat[[x]][sort.int(log2mat1RegionRowMeans,
                           decreasing = T,
                           index.return = T,
                           na.last = T)$ix,]
}, mc.cores = length(log2covMat))

# Convert matrices to class "normalizedMatrix" with associated attributes
# for use by EnrichedHeatmap function
for(x in seq_along(log2covMatSorted)) {
  attr(log2covMatSorted[[x]], "upstream_index") = 1:(upstream/binSize)
  attr(log2covMatSorted[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(log2covMatSorted[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(log2covMatSorted[[x]], "extend") = c(upstream, downstream)
  attr(log2covMatSorted[[x]], "smooth") = TRUE
  attr(log2covMatSorted[[x]], "signal_name") = libNamePlot[x]
  attr(log2covMatSorted[[x]], "target_name") = featureName
  attr(log2covMatSorted[[x]], "target_is_single_point") = FALSE
  attr(log2covMatSorted[[x]], "background") = 0
  attr(log2covMatSorted[[x]], "signal_is_categorical") = FALSE
  class(log2covMatSorted[[x]]) = c("normalizedMatrix", "matrix")
}

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Heatmap plotting function
featureHeatmap <- function(matSorted,
                           col_fun,
                           colour,
                           datName,
                           rowOrder) {
  print(EnrichedHeatmap(mat = matSorted,
                        col = col_fun,
                        top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = colour,
                                                                                              lwd = 3),
                                                                                    yaxis_side= "right",
                                                                                    yaxis_facing = "right",
                                                                                    yaxis_gp =gpar(fontsize = 10),
                                                                                    pos_line_gp = gpar(col = "black",
                                                                                                       lty = 2,
                                                                                                       lwd = 2))),
                        top_annotation_height = unit(2, "cm"),
                        name = "Coverage",
                        row_order = rowOrder,
                        column_title = datName,
                        axis_name = c(paste0("-", flankNamePlot),
                                      featureStartLab, featureEndLab,
                                      paste0("+", flankNamePlot)),
                        axis_name_gp = gpar(fontsize = 12),
                        border = FALSE,
                        pos_line_gp = gpar(col = "black", lty = 2, lwd = 2),
                        use_raster = FALSE))
}

# Define heatmap colours
rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")

# Plot
# ChIP-seq
foreach(x = 1:length(libNameChIP)) %dopar% {
  log2ChIP_col_fun <- colorRamp2(quantile(c(log2covMatSorted[[x]]),
                                          c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95),
                                          na.rm = T),
                                 rich8to6equal)
  pdf(paste0(plotDir, "log2_", libNameChIP[x], "_control_around_", featureName,
             "_heatmap_ordered_by_", libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  featureHeatmap(matSorted = log2covMatSorted[[x]],
                 col_fun = log2ChIP_col_fun,
                 colour = "red",
                 datName = libNamePlot[x],
                 rowOrder = c(1:dim(log2covMatSorted[[x]])[1]))
  dev.off()
}
