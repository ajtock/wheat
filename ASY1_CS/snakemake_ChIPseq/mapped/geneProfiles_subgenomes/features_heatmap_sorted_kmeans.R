#!/applications/R/R-3.4.0/bin/Rscript

# Cluster genes by log2(ChIP/control)

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_heatmap_sorted_kmeans.R ASY1_CS_Rep1_ChIP ASY1_CS both genes_in_Agenome_genomewide 3500 2000 2kb '2 kb' 20 20bp terminators

libName <- "ASY1_CS_Rep1_ChIP"
dirName <- "ASY1_CS"
align <- "both"
featureName <- "genes_in_Agenome_genomewide"
bodyLength <- 3500
upstream <- 2000
downstream <- 2000
flankName <- "2kb"
flankNamePlot <- "2 kb"
binSize <- 20
binName <- "20bp"
region <- "terminators"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
align <- args[3]
featureNumber <- as.numeric(args[4])
featureName <- args[5]
bodyLength <- as.numeric(args[6])
upstream <- as.numeric(args[7])
downstream <- as.numeric(args[7])
flankName <- args[8]
flankNamePlot <- args[9]
binSize <- as.numeric(args[10])
binName <- args[11]
region <- args[12]

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
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

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
controlDirs <- c(
                 "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/input/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                 "/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/"
                )
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

# Determine number of clusters to be used for k-means clustering
# by generating a scree ("elbow") plot of the ratio of the
# within-cluster sum of squares (WSS) to the total sum of squares (TSS)
# for k clusters

# First run garbage collection to free up memory
gc()
# Set seed for reproducible results
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
pdf(paste0(outDir, "screePlot_k_clusters_log2_",
           libName, "_control_around_",
           region, "_of_", featureName, ".pdf"))
plot(x = 1:kMax, y = ratio_ss,
     type = "b",
     xlab = "k", ylab = "WSS / TSS")
abline(h = 0.2, lty = 2, col = "red")
dev.off()

# First run garbage collection to free up memory
gc()
# Set seed for reproducible results
set.seed(938402845)

# 4 clusters
kDef <- 4
km <- kmeans(x = log2ChIPmatRegion,
             centers = kDef,
             iter.max = 10,
             nstart = 10)

## Calculate Dunn's index to assess compactness and separation of clusters
## "The Dunn Index is the ratio of the smallest distance between observations
## not in the same cluster to the largest intra-cluster distance.
## The Dunn Index has a value between zero and infinity, and should be maximized."
## (NOTE: takes a long time)
#library(clValid)
#km_dunn <- dunn(clusters = km$cluster, Data = log2ChIPmatRegion)

# Get feature indices for each cluster
featureIndicesList <- lapply(seq_along(1:kDef), function(k) {
  which(km$cluster == k)
})

# Load features 
features <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                              substring(featureName, first = 10), ".gff3"), header = F)
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
                            "_by_log2_", libName, "_control_in_", region, ".txt"),
                            quote = F, row.names = F, col.names = F)
})


