#!/applications/R/R-3.5.0/bin/Rscript

# Plot boxplots and violin plots of recombination rate for each gene quantile; e.g.,
# quantiles_by_log2_DMC1_Rep1_ChIP_control_in_promoters/features_4quantiles_by_log2_DMC1_Rep1_ChIP_control_in_promoters_of_genes_in_Agenome_Bgenome_Dgenome_genomewide.tx

# Usage:
# /applications/R/R-3.5.0/bin/Rscript quantile_genes_cMMb_density_plot.R DMC1_Rep1_ChIP DMC1 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' promoters 4

#libName <- "DMC1_Rep1_ChIP"
#dirName <- "DMC1"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#region <- "promoters"
#quantiles <- 4

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
region <- args[4]
quantiles <- as.numeric(args[5])

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
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
featureNamePlot <- paste0(sub("_\\w+", "", dirName), " ",
                          substr(featureName[1], start = 1, stop = 4),
                          " quantiles")
ranFeatNamePlot <- paste0("Random ",
                          substr(featureName[1], start = 1, stop = 4),
                          " quantiles")
ranLocNamePlot <- "Random locus quantiles"

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")
makeTransparent <- function(thisColour, alpha = 180)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
quantileColours <- makeTransparent(quantileColours)

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
# Load table of ranLocs grouped according to feature quantiles
ranLocsDF <- read.table(paste0(outDir, "features_", quantiles, "quantiles",
                               "_by_log2_", libName, "_control_in_",
                               region, "_of_",
                               substring(featureName[1][1], first = 1, last = 5), "_in_",
                               paste0(substring(featureName, first = 10, last = 16),
                                      collapse = "_"), "_",
                               substring(featureName[1][1], first = 18), "_ranLocs.txt"),
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
lapply(seq_along(1:quantiles), function(k) {
  sapply(seq_along(chrs), function(x) {
    if(!identical(dim(featuresDF[randomPCIndices[[k]],][featuresDF[randomPCIndices[[k]],]$seqnames == chrs[x],]),
                  dim(featuresDF[quantileIndices[[k]],][featuresDF[quantileIndices[[k]],]$seqnames == chrs[x],])))     {
      stop("Quantile features and random features do not consist of the same number of features per chromosome")
    }
  })
})

featuresDFtmp <- data.frame(featuresDF,
                            random = as.character(""),
                            stringsAsFactors = F)
ranFeatsDF <- data.frame()
for(k in 1:quantiles) {
  featuresDFtmp[randomPCIndices[[k]],]$random <- paste0("Random ", k)
  ranFeatsDFk <- featuresDFtmp[featuresDFtmp$random == paste0("Random ", k),]
  ranFeatsDF <- rbind(ranFeatsDF, ranFeatsDFk)
}

xmin <- min(c(
              featuresDF[unlist(quantileIndices),]$cMMb,
              featuresDF[unlist(randomPCIndices),]$cMMb,
              ranLocsDF[unlist(quantileIndices),]$cMMb
             ), na.rm = T)
xmax <- max(c(
              featuresDF[unlist(quantileIndices),]$cMMb,
              featuresDF[unlist(randomPCIndices),]$cMMb,
              ranLocsDF[unlist(quantileIndices),]$cMMb
             ), na.rm = T)
xmin <- 0
xmax <- 3
minDensity <- 0
maxDensity <- max(density(featuresDF[featuresDF$quantile == "Quantile 4",]$cMMb,
                          na.rm = T)$y)+2
maxDensity <- max(
  c(
    sapply(1:quantiles, function(k) {
#      max(c(max(density(ranLocsDF[ranLocsDF$random == paste0("Random ", k),]$cMMb,
#                        na.rm = T)$y),
      max(c(max(density(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$cMMb,
                        na.rm = T)$y),
            max(density(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),]$cMMb,
                        na.rm = T)$y)))
     })
   )
)+2

# Recombination rate (cM/Mb) boxplot or violin plot function
cMMb_plotFun <- function(featuresDF,
                         parameter,
                         parameterLab,
                         featureGroup,
                         featureNamePlot,
                         quantileColours) {
  ggplot(data = featuresDF,
         mapping = aes(x = get(parameter),
                       colour = reorder(x = get(featureGroup), X = desc(get(featureGroup))),
                       group = reorder(x = get(featureGroup), X = desc(get(featureGroup))))) +
  scale_colour_manual(values = rev(quantileColours)) +
  geom_density(size = 1.5) +
  scale_x_continuous(limits = c(xmin, xmax),
                     labels = function(x) sprintf("%1.1f", x)) +
  scale_y_continuous(limits = c(minDensity, maxDensity),
                     labels = function(x) sprintf("%1.1f", x)) +
  labs(x = parameterLab,
       y = "Density") +
  theme_bw() +
  theme(axis.line.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.text.x = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 26, colour = "black"),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 22, colour = "black"),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.1,0.3),"cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot)))
}

ggObjGA_feature <- cMMb_plotFun(featuresDF = featuresDF,
                                parameter = "cMMb",
                                parameterLab = "Recombination rate (cM/Mb)",
                                featureGroup = "quantile", 
                                featureNamePlot = featureNamePlot,
                                quantileColours = quantileColours
                               )
ggObjGA_ranFeat <- cMMb_plotFun(featuresDF = ranFeatsDF,
                                parameter = "cMMb",
                                parameterLab = "Recombination rate (cM/Mb)",
                                featureGroup = "random", 
                                featureNamePlot = ranFeatNamePlot,
                                quantileColours = quantileColours
                               )
ggObjGA_ranLocs <- cMMb_plotFun(featuresDF = ranLocsDF,
                                parameter = "cMMb",
                                parameterLab = "Recombination rate (cM/Mb)",
                                featureGroup = "random", 
                                featureNamePlot = ranLocNamePlot,
                                quantileColours = quantileColours
                               )
ggObjGA_combined <- grid.arrange(ggObjGA_feature,
                                 ggObjGA_ranFeat,
#                                 ggObjGA_ranLocs,
                                 ncol = 2, as.table = F)
ggsave(paste0(plotDir,
              "cMMb_around_", quantiles, "quantiles",
              "_by_log2_", libName, "_control_in_", region, "_of_",
              substring(featureName[1][1], first = 1, last = 5), "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v231119.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 14)

