#!/applications/R/R-3.4.0/bin/Rscript

# Plot average coverage profiles with 95% CIs around
# Cluster 1 genes annotated with the  "defense response" GO term; i.e.,
# clusters_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters/GO/cluster1_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters_of_genes_in_Agenome_genomewide_GO_BP/cluster1_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters_of_genes_in_Agenome_genomewide_GO_BP_enrichment_GO\:0006952.txt

# Usage:
# /applications/R/R-3.4.0/bin/Rscript defense_response_genes_avgProfileRibbon.R genes_in_Agenome_genomewide 3500 2000 2kb '2 kb' 20 20bp promoters '1' '4' '0006952' both ASY1_CS_Rep1_ChIP ASY1_CS purple4

featureName <- "genes_in_Agenome_genomewide"
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
libName <- "ASY1_CS_Rep1_ChIP"
dirName <- "ASY1_CS"
libNamePlot <- "ASY1"
colour <- "purple4"

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
bodyLength <- as.numeric(args[2])
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[4])
flankName <- args[5]
flankNamePlot <- args[6]
binSize <- as.numeric(args[7])
binName <- args[8]
region <- args[9]
clusterNo <- as.character(args[10])
clusterLast <- as.character(args[11])
GO_ID <- as.character(args[12])
align <- as.character(args[13])
libName <- args[14]
dirName <- args[15]
libNamePlot <- args[16]
colour <- args[17]

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

# Add column names
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

## feature
# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  data.frame(window = colnames(log2ChIP_featureMats[[x]]),
             t(log2ChIP_featureMats[[x]]))
}, mc.cores = length(log2ChIP_featureMats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  gather(data  = wideDFfeature_list_log2ChIP[[x]],
         key   = feature,
         value = coverage,
         -window)
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  tidyDFfeature_list_log2ChIP[[x]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]]$window,
                                                    levels = as.character(wideDFfeature_list_log2ChIP[[x]]$window))
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

# Convert list summaryDFfeature_list_log2ChIP into a single data.frame for plotting
summaryDFfeature <- bind_rows(summaryDFfeature_list_log2ChIP, .id = "libName")
summaryDFfeature$libName <- factor(summaryDFfeature$libName,
                                   levels = names(summaryDFfeature_list_log2ChIP))


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
