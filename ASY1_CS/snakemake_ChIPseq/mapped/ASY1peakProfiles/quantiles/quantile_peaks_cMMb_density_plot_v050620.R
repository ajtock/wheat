#!/applications/R/R-3.5.0/bin/Rscript

# Plot density and means with 95% CIs of recombination rate for each peak quantile

# Usage:
# /applications/R/R-3.5.0/bin/Rscript quantile_peaks_cMMb_density_plot_v050620.R cMMb ASY1_CS_Rep1_ChIP ASY1_CS 'Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin' 4

#orderingFactor <- "cMMb"
#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#featureName <- unlist(strsplit("Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin",
#                               split = ","))
#quantiles <- 4

args <- commandArgs(trailingOnly = T)
orderingFactor <- args[1]
libName <- args[2]
dirName <- args[3]
featureName <- unlist(strsplit(args[4],
                               split = ","))
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

if(orderingFactor == "cMMb") {
outDir <- paste0("quantiles_by_", orderingFactor, "/")
} else {
outDir <- paste0("quantiles_by_log2_", orderingFactor,
                 "_control_in_", region, "/")
}
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Define plot titles
if(orderingFactor == "cMMb") {
  featureNamePlot <- paste0("cM/Mb ",
                            gsub("_\\w+$", "", libName),
                            " peak quantiles")
}
ranFeatNamePlot <- paste0("Random ",
                          gsub("_\\w+$", "", libName),
                          " peak quantiles")
#ranLocNamePlot <- "Random locus quantiles"

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")
makeTransparent <- function(thisColour, alpha = 250)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
quantileColours <- makeTransparent(quantileColours)

# Genomic definitions
#chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
#chrs <- chrs[-length(chrs)]
chrs <- paste0(rep("chr", 21), rep(1:7, 3),
               c(rep("A", 7), rep("B", 7), rep("D", 7)))

# Load table of features grouped into quantiles
# by decreasing cM/Mb
featuresDF <- read.table(paste0(outDir,
                                "features_", quantiles, "quantiles",
                                 "_by_", orderingFactor,
                                 "_of_", libName, "_peaks_in_",
                                 paste0(featureName,
                                        collapse = "_"), ".txt"),
                         header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
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
# If features from multiple subgenomes and/or compartments are to be analysed,
# concatenate the corresponding feature data.frames
if(length(featureName) > 1) {
  features <- do.call(rbind, features)
} else {
  features <- features[[1]]
}
colnames(features) <- c("chr", "start", "end", "name", "score", "strand", "featureID")
stopifnot(identical(as.character(featuresDF$featureID),
                    as.character(features$featureID)))
rm(features); gc()

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
#set.seed(453838430)


# Divide features into random sets of equal number,
# with the same number of genes per chromosome as
# above-defined orderingFactor-defined feature quantiles
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

# Calculate means, SDs, SEMs and 95% CIs
# and create dataframe of summary statistics for plotting
featuresDF_quantileMean <- sapply(1:quantiles, function(k) {
  mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$cMMb, na.rm = T)
})
featuresDF_quantileSD <- sapply(1:quantiles, function(k) {
  sd(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$cMMb, na.rm = T)
})
featuresDF_quantileSEM <- sapply(1:quantiles, function(k) {
  featuresDF_quantileSD[k] / sqrt( (dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1] - 1) )
})
featuresDF_quantileCIlower <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] -
    ( qt(0.975, df = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_quantileCIupper <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] +
    ( qt(0.975, df = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_summary_stats <- data.frame(quantile = paste0("Quantile ", 1:4),
                                       Mean = featuresDF_quantileMean,
                                       SD = featuresDF_quantileSD,
                                       SEM = featuresDF_quantileSEM,
                                       CIlower = featuresDF_quantileCIlower,
                                       CIupper = featuresDF_quantileCIupper,
                                       stringsAsFactors = F)

ranFeatsDF_randomMean <- sapply(1:quantiles, function(k) {
  mean(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),]$cMMb, na.rm = T)
})
ranFeatsDF_randomSD <- sapply(1:quantiles, function(k) {
  sd(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),]$cMMb, na.rm = T)
})
ranFeatsDF_randomSEM <- sapply(1:quantiles, function(k) {
  ranFeatsDF_randomSD[k] / sqrt( (dim(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),])[1] - 1) )
})
ranFeatsDF_randomCIlower <- sapply(1:quantiles, function(k) {
  ranFeatsDF_randomMean[k] -
    ( qt(0.975, df = dim(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),])[1]-1 ) *
      ranFeatsDF_randomSEM[k] )
})
ranFeatsDF_randomCIupper <- sapply(1:quantiles, function(k) {
  ranFeatsDF_randomMean[k] +
    ( qt(0.975, df = dim(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),])[1]-1 ) *
      ranFeatsDF_randomSEM[k] )
})
ranFeatsDF_summary_stats <- data.frame(random = paste0("Random ", 1:4),
                                       Mean = ranFeatsDF_randomMean,
                                       SD = ranFeatsDF_randomSD,
                                       SEM = ranFeatsDF_randomSEM,
                                       CIlower = ranFeatsDF_randomCIlower,
                                       CIupper = ranFeatsDF_randomCIupper,
                                       stringsAsFactors = F)
summary_stats_min <- min(c(featuresDF_summary_stats$CIlower, ranFeatsDF_summary_stats$CIlower), na.rm = T)
summary_stats_max <- max(c(featuresDF_summary_stats$CIupper, ranFeatsDF_summary_stats$CIupper), na.rm = T)

# Take top 95% of data to aid visualisation in density plots
featuresDF <- featuresDF[which(featuresDF$cMMb <=
                               quantile(featuresDF$cMMb,
                                        probs = 0.99, na.rm = T)),]
ranFeatsDF <- ranFeatsDF[which(ranFeatsDF$cMMb <=
                               quantile(ranFeatsDF$cMMb,
                                        probs = 0.99, na.rm = T)),]

xmin <- min(c(
              featuresDF[unlist(quantileIndices),]$cMMb,
              featuresDF[unlist(randomPCIndices),]$cMMb
             ), na.rm = T)
xmax <- max(c(
              featuresDF[unlist(quantileIndices),]$cMMb,
              featuresDF[unlist(randomPCIndices),]$cMMb
             ), na.rm = T)
minDensity <- 0
maxDensity <- max(density(featuresDF[featuresDF$quantile == "Quantile 4",]$cMMb,
                          na.rm = T)$y)+0.04
maxDensity <- max(
  c(
    sapply(1:quantiles, function(k) {
      max(c(max(density(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$cMMb,
                        na.rm = T)$y),
            max(density(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),]$cMMb,
                        na.rm = T)$y)))
     })
   )
)+0.04

# Define legend labels
legendLabs_feature <- lapply(1:quantiles, function(x) {
  grobTree(textGrob(bquote(.(paste0("Quantile ", 1:quantiles)[x])),
                    x = 0.65, y = 0.90-((x-1)*0.07), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 22)))
})
legendLabs_ranFeat <- lapply(1:quantiles, function(x) {
  grobTree(textGrob(bquote(.(paste0("Random ", 1:quantiles)[x])),
                    x = 0.65, y = 0.90-((x-1)*0.07), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 22)))
})

# Recombination rate (cM/Mb) density plot function
cMMb_plotFun <- function(lociDF,
                         parameter,
                         parameterLab,
                         featureGroup,
                         featureNamePlot,
                         legendLabs,
                         quantileColours) {
  ggplot(data = lociDF,
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
  annotation_custom(legendLabs[[1]]) +
  annotation_custom(legendLabs[[2]]) +
  annotation_custom(legendLabs[[3]]) +
  annotation_custom(legendLabs[[4]]) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 2.0, colour = "black"),
        axis.line.x = element_line(size = 2.0, colour = "black"),
        axis.ticks.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.x = element_line(size = 2.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.text.x = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 26, colour = "black"),
        legend.position = "none",
##       TEST WITH THIS WAY OF PLOTTING LEGEND LABELS TO CHECK THAT QUANTILES CORRESPOND TO EXPECTED COLOURS
#        legend.position = c(0.8, 0.8),
#        legend.text = element_text(size = 22, colour = "black"),
#        legend.key.size = unit(1, "cm"),
#        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.3,0.3),"cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot)))
}

# Plot means and 95% confidence intervals
cMMb_meanCIs <- function(dataFrame,
                         parameterLab,
                         featureGroup,
                         featureNamePlot,
                         quantileColours) {
  ggplot(data = dataFrame,
         mapping = aes(x = get(featureGroup),
                       y = Mean,
                       colour = get(featureGroup))) +
  labs(colour = "") +
  geom_point(shape = 19, size = 6, position = position_dodge(width = 0.2)) +
  geom_errorbar(mapping = aes(ymin = CIlower,
                              ymax = CIupper),
                width = 0.2, size = 2, position = position_dodge(width = 0.2)) +
  scale_colour_manual(values = quantileColours) +
  scale_y_continuous(limits = c(summary_stats_min, summary_stats_max),
                     labels = function(x) sprintf("%1.2f", x)) +
#  scale_x_discrete(breaks = as.vector(dataFrame$quantile),
#                   labels = as.vector(dataFrame$quantile)) +
  labs(x = "",
       y = parameterLab) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.text.x = element_text(size = 22, colour = quantileColours, hjust = 1.0, vjust = 1.0, angle = 45),
        axis.title = element_text(size = 26, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.1,0.3),"cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot)))
}

ggObjGA_feature <- cMMb_plotFun(lociDF = featuresDF,
                                parameter = "cMMb",
                                parameterLab = "Recombination rate (cM/Mb)",
                                featureGroup = "quantile", 
                                featureNamePlot = featureNamePlot,
                                legendLabs = legendLabs_feature,
                                quantileColours = quantileColours
                               )
ggObjGA_ranFeat <- cMMb_plotFun(lociDF = ranFeatsDF,
                                parameter = "cMMb",
                                parameterLab = "Recombination rate (cM/Mb)",
                                featureGroup = "random", 
                                featureNamePlot = ranFeatNamePlot,
                                legendLabs = legendLabs_ranFeat,
                                quantileColours = quantileColours
                               )
ggObjGA_feature_mean <- cMMb_meanCIs(dataFrame = featuresDF_summary_stats,
                                     parameterLab = "Recombination rate (cM/Mb)",
                                     featureGroup = "quantile",
                                     featureNamePlot = featureNamePlot,
                                     quantileColours = quantileColours
                                    )
ggObjGA_ranFeat_mean <- cMMb_meanCIs(dataFrame = ranFeatsDF_summary_stats,
                                     parameterLab = "Recombination rate (cM/Mb)",
                                     featureGroup = "random",
                                     featureNamePlot = ranFeatNamePlot,
                                     quantileColours = quantileColours
                                    )
ggObjGA_combined <- grid.arrange(ggObjGA_feature,
                                 ggObjGA_feature_mean,
                                 ggObjGA_ranFeat,
                                 ggObjGA_ranFeat_mean,
                                 ncol = 2, as.table = F)
ggsave(paste0(plotDir,
              "cMMb_around_", quantiles, "quantiles",
              "_by_", orderingFactor,
              "_of_", libName, "_peaks_in_",
              paste0(featureName,
                     collapse = "_"),
              "_v050620.pdf"),
       plot = ggObjGA_combined,
       height = 13, width = 14)
