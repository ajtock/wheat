#!/applications/R/R-3.4.0/bin/Rscript

# Plot heatmaps of features sorted by coverage levels in bodies,
# within feature bodies, or in terminators

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_heatmap_sorted.R ASY1_CS_Rep1_ChIP ASY1_CS both ASY1_CS_peaks_in_Agenome_euchromatin 400 2000 2kb '2 kb' 20 20bp bodies

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#align <- "both"
#featureName <- "ASY1_CS_peaks_in_Agenome_euchromatin"
#bodyLength <- 400
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "bodies"

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

plotDir <- "./plots/"
regionPlotDir <- paste0(plotDir, "log2_", libName, "_control_in_", region, "/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
system(paste0("[ -d ", regionPlotDir, " ] || mkdir ", regionPlotDir))
plotDir <- regionPlotDir

# Load libName matrix, calculate log2(ChIP/control), and extract region for sorting of features
mat1 <- as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/",
                                    dirName,
                                    "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
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
           controlNamesDir[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else if(controlNames[x] == "MNase_Rep1") {
    paste0("/home/ajt200/analysis/wheat/",
           controlNamesDir[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
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
# for each matrix depending on library
log2mat1 <- if(libName %in% c(
                              "ASY1_CS_Rep1_ChIP",
                              "DMC1_Rep1_ChIP",
                              "H3K4me3_ChIP_SRR6350668",
                              "H3K27me3_ChIP_SRR6350666",
                              "H3K36me3_ChIP_SRR6350670",
                              "H3K9ac_ChIP_SRR6350667",
                              "H3K4me1_Rep1_ChIP_SRR8126618",
                              "H3K27ac_Rep1_ChIP_SRR8126621"
                             )) {
  print(paste0(libName, " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
  log2((mat1+1)/(controlmats[[1]]+1))
} else {
  print(paste0(libName, " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
  log2((mat1+1)/(controlmats[[2]]+1))
}

# Extract region for sorting of features by decreasing log2mat1 levels (adjust promoter/terminator size as necessary)
if( region == "promoters" ) {
  log2mat1Region <- log2mat1[,(((upstream-500)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  log2mat1Region <- log2mat1[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(500/binSize))]
} else if ( region == "bodies" ) {
  log2mat1Region <- log2mat1[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else {
  print("The region name provided does not match 'promoters', 'terminators', or 'bodies'")
}
log2mat1RegionRowMeans <- rowMeans(log2mat1Region, na.rm = T)
log2mat1RegionRowMeansSorted <- sort.int(log2mat1RegionRowMeans,
                                         decreasing = T,
                                         index.return = T,
                                         na.last = T)
log2mat1RegionRowSums <- rowSums(log2mat1Region, na.rm = T)
log2mat1RegionRowSumsSorted <- sort.int(log2mat1RegionRowSums,
                                        decreasing = T,
                                        index.return = T,
                                        na.last = T)
log2mat1RegionRowMedians <- apply(X = log2mat1Region,
                                  MARGIN = 1,
                                  FUN = median)
log2mat1RegionRowMediansSorted <- sort.int(log2mat1RegionRowMedians,
                                           decreasing = T,
                                           index.return = T,
                                           na.last = T)
log2mat1Sorted <- log2mat1[sort.int(log2mat1RegionRowMeans,
                                    decreasing = T,
                                    index.return = T,
                                    na.last = T)$ix,]

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
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else if(ChIPNames[x] %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                                "H3K27ac_Rep1_ChIP_SRR8126621")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else {
    paste0("/home/ajt200/analysis/wheat/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  }
})
#ChIPDirs <- c(
#              "/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
#              "/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
#              "/home/ajt200/analysis/wheat/H2AZ/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
#              "/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
#              "/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/H3K4me1/snakemake_ChIPseq/mapped/matrices/",
#              "/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/H3K27ac/snakemake_ChIPseq/mapped/matrices/",
#              "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
#              "/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
#              "/home/ajt200/analysis/wheat/H3K27me1/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/"
#             )
otherDirs <- sapply(seq_along(otherNames), function(x) {
  if(otherNames[x] %in% c("MNase_Rep1")) {
    paste0("/home/ajt200/analysis/wheat/",
           otherNamesDir[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else if(otherNames[x] %in% c("DNaseI_Rep1_SRR8447247")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
           otherNamesDir[x], "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
  } else {
    if(!(otherNames %in% c("MNase_Rep1", "DNaseI_Rep1_SRR8447247"))) {
      stop("otherNames[x] is neither MNase_Rep1 nor DNaseI_Rep1_SRR8447247")
    }
  }
})
#otherDirs <- c(
#               "/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/",
#               "/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/DNase1/snakemake_ChIPseq/mapped/matrices/"
#              )
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

# Sorted by decreasing log2mat1RegionRowMeans (e.g., ASY1 in feature bodies)
log2ChIPmatsSorted <- mclapply(seq_along(log2ChIPmats), function(x) {
  log2ChIPmats[[x]][sort.int(log2mat1RegionRowMeans,
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

othermatsSorted <- mclapply(seq_along(othermats), function(x) {
  othermats[[x]][sort.int(log2mat1RegionRowMeans,
                          decreasing = T,
                          index.return = T,
                          na.last = T)$ix,]
}, mc.cores = length(othermats))

for(x in seq_along(othermatsSorted)) {
  attr(othermatsSorted[[x]], "upstream_index") = 1:(upstream/binSize)
  attr(othermatsSorted[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(othermatsSorted[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(othermatsSorted[[x]], "extend") = c(upstream, downstream)
  attr(othermatsSorted[[x]], "smooth") = FALSE
  attr(othermatsSorted[[x]], "signal_name") = ChIPNamesPlot[x]
  attr(othermatsSorted[[x]], "target_name") = featureName
  attr(othermatsSorted[[x]], "target_is_single_point") = FALSE
  attr(othermatsSorted[[x]], "background") = 0
  attr(othermatsSorted[[x]], "signal_is_categorical") = FALSE
  class(othermatsSorted[[x]]) = c("normalizedMatrix", "matrix")
}

controlmatsSorted <- mclapply(seq_along(controlmats), function(x) {
  controlmats[[x]][sort.int(log2mat1RegionRowMeans,
                            decreasing = T,
                            index.return = T,
                            na.last = T)$ix,]
}, mc.cores = length(controlmats))

for(x in seq_along(controlmatsSorted)) {
  attr(controlmatsSorted[[x]], "upstream_index") = 1:(upstream/binSize)
  attr(controlmatsSorted[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(controlmatsSorted[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(controlmatsSorted[[x]], "extend") = c(upstream, downstream)
  attr(controlmatsSorted[[x]], "smooth") = FALSE
  attr(controlmatsSorted[[x]], "signal_name") = ChIPNamesPlot[x]
  attr(controlmatsSorted[[x]], "target_name") = featureName
  attr(controlmatsSorted[[x]], "target_is_single_point") = FALSE
  attr(controlmatsSorted[[x]], "background") = 0
  attr(controlmatsSorted[[x]], "signal_is_categorical") = FALSE
  class(controlmatsSorted[[x]]) = c("normalizedMatrix", "matrix")
}

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Heatmap plotting function
# Note that for plotting heatmaps for individual datasets in separate PDFs,
# must edit this function - print(EnrichedHeatmap(...))
featureHeatmap <- function(matSorted,
                           col_fun,
                           colour,
                           datName,
                           rowOrder) {
  EnrichedHeatmap(mat = matSorted,
                  col = col_fun,
                  row_order = rowOrder,
                  column_title = datName,
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = colour,
                                                                                        lwd = 3),
                                                                              yaxis_side = "right",
                                                                              yaxis_facing = "right",
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
rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")

# Plot together
# ChIP
log2ChIPhtmpList <- mclapply(seq_along(ChIPNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(log2ChIPmatsSorted[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(matSorted = log2ChIPmatsSorted[[x]],
                 col_fun = ChIP_col_fun,
                 colour = ChIPColours[x],
                 datName = ChIPNamesPlot[x],
                 rowOrder = c(1:dim(log2ChIPmatsSorted[[x]])[1]))
}, mc.cores = length(log2ChIPmatsSorted))
otherhtmpList <- mclapply(seq_along(otherNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(othermatsSorted[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(matSorted = othermatsSorted[[x]],
                 col_fun = ChIP_col_fun,
                 colour = otherColours[x],
                 datName = otherNamesPlot[x],
                 rowOrder = c(1:dim(othermatsSorted[[x]])[1]))
}, mc.cores = length(othermatsSorted))
controlhtmpList <- mclapply(seq_along(controlNames), function(x) {
  ChIP_col_fun <- colorRamp2(quantile(controlmatsSorted[[x]],
                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                      na.rm = T),
                             rich8to6equal)
  featureHeatmap(matSorted = controlmatsSorted[[x]],
                 col_fun = ChIP_col_fun,
                 colour = controlColours[x],
                 datName = controlNamesPlot[x],
                 rowOrder = c(1:dim(controlmatsSorted[[x]])[1]))
}, mc.cores = length(controlmatsSorted))

htmpList <- c(log2ChIPhtmpList, otherhtmpList, controlhtmpList[[1]])

htmps <- NULL
for(x in 1:length(htmpList)) {
  htmps <- htmps + htmpList[[x]]
}
pdf(paste0(plotDir, "log2ChIPcontrol_around_", featureName,
           "_heatmaps_ordered_by_", libName, "_in_", region, ".pdf"),
    width = 3*length(htmpList),
    height = 8)
draw(htmps,
     heatmap_legend_side = "bottom",
     gap = unit(c(14), "mm"))
dev.off()
