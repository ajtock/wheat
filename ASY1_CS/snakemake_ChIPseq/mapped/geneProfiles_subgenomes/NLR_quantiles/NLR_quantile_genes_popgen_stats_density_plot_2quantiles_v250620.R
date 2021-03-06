#!/applications/R/R-3.5.0/bin/Rscript

# Plot density profiles and means and 95% CIs of population genetics statistics for each gene quantile

# Usage:
# /applications/R/R-3.5.0/bin/Rscript NLR_quantile_genes_popgen_stats_density_plot_2quantiles_v250620.R HudsonRM_all 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' genes 2 TajimaD_all "Tajima's D" 1.0 0.2 '%1.1f' '%1.1f' '%1.1f' 0.65

# /applications/R/R-3.5.0/bin/Rscript NLR_quantile_genes_popgen_stats_density_plot_2quantiles_v250620.R HudsonRM_all 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' genes 2 RozasR2_all "Rozas' R 2" 0.95 0.2 '%1.2f' '%2.0f' '%2.0f' 0.38

# /applications/R/R-3.5.0/bin/Rscript NLR_quantile_genes_popgen_stats_density_plot_2quantiles_v250620.R HudsonRM_all 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' genes 2 CLR_all "CLR" 1.0 0.005 '%1.0f' '%1.2f' '%1.0f' 0.65

# /applications/R/R-3.5.0/bin/Rscript NLR_quantile_genes_popgen_stats_density_plot_2quantiles_v250620.R HudsonRM_all 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' genes 2 HudsonRM_all "Hudson's R M" 0.90 0.2 '%1.2f' '%2.0f' '%1.2f' 0.38

# /applications/R/R-3.5.0/bin/Rscript NLR_quantile_genes_popgen_stats_density_plot_2quantiles_v250620.R HudsonRM_all 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' genes 2 nucleotideDiversity_all "Diversity pi" 0.95 0.2 '%1.2f' '%3.0f' '%1.2f' 0.65

# /applications/R/R-3.5.0/bin/Rscript NLR_quantile_genes_popgen_stats_density_plot_2quantiles_v250620.R HudsonRM_all 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' genes 2 nSegregatingSites_all "Norm. seg. sites" 1.0 0.05 '%1.1f' '%1.1f' '%1.1f' 0.65

#libName <- "HudsonRM_all"
#featureName <- unlist(strsplit("NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide",
#                               split = ","))
#region <- "genes"
#quantiles <- 2
#orderingFactor <- "TajimaD_all"
##orderingFactor <- "RozasR2_all"
#orderingFactorName <- bquote("Tajima's" ~ italic("D"))
##orderingFactorName <- bquote("Rozas'" ~ italic("R")[2])
#orderingFactorName <- unlist(strsplit("Tajima's D", split = " "))
##orderingFactorName <- unlist(strsplit("Rozas' R 2", split = " "))
#densityProp <- 0.99
#maxDensityPlus <- 0.2
#xDec <- "%1.1f"
#yDec <- "%1.1f"
#yDec2 <- "%1.2f"
#legendLabX <- "0.65"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
featureName <- unlist(strsplit(args[2],
                               split = ","))
region <- args[3]
quantiles <- as.numeric(args[4])
orderingFactor <- args[5]
orderingFactorName <- unlist(strsplit(args[6], split = " "))
if(grepl("Tajima", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ italic(.(orderingFactorName[2])))
} else if(grepl("Rozas' Z", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ italic(.(orderingFactorName[2])))
} else if(grepl("Rozas' R", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ italic(.(orderingFactorName[2]))[.(as.numeric(orderingFactorName[3]))])
} else if(grepl("Hudson's R", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ italic(.(orderingFactorName[2])[.(as.character(orderingFactorName[3]))]))
} else if(grepl("Diversity", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ "(" * .(as.symbol(orderingFactorName[2])) * ")")
} else {
  orderingFactorName <- paste(orderingFactorName, collapse = " ")
}
densityProp <- as.numeric(args[7])
maxDensityPlus <- as.numeric(args[8])
xDec <- as.character(args[9])
yDec <- as.character(args[10])
yDec2 <- as.character(args[11])
legendLabX <- as.numeric(args[12])

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(ggtext)
library(grid)
library(gridExtra)
library(extrafont)

pop_name_plot <- c(
                   "Africa",
                   "Middle East",
                   "Asia",
                   "Former SU",
                   "Eastern Europe",
                   "Western Europe",
                   "North America",
                   "Central America",
                   "South America",
                   "Oceania"
                  )
pop_name <- gsub(" ", "", pop_name_plot)

if(libName %in% c("cMMb", "cluster_members", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  outDir <- paste0("quantiles_by_", libName, "/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", region, "/")
}
outDir <- sapply(seq_along(pop_name), function(x) {
  paste0(outDir, pop_name[x], "/")
})
plotDir <- paste0(outDir, "plots/")

# Define plot titles
if(libName %in% "cMMb") {
  featureNamePlot <- paste0("cM/Mb ",
                            substr(featureName[1], start = 1, stop = 3),
                            " quantiles")
} else if(libName %in% "cluster_members") {
  featureNamePlot <- paste0(
                            substr(featureName[1], start = 1, stop = 3),
                            " cluster size",
                            " quantiles")
} else if(libName %in% c("HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  featureNamePlot <- bquote("Hudson's" ~ italic("R"[M]) ~
                            .(substr(featureName[1], start = 1, stop = 3)) *
                            " quantiles")
} else {
  featureNamePlot <- paste0(sub("_\\w+", "", libName), " ",
                            substr(featureName[1], start = 1, stop = 3),
                            " quantiles (", region, ")")
}
ranFeatNamePlot <- paste0("Random ",
                          substr(featureName[1], start = 1, stop = 3),
                          " quantiles")
#ranLocNamePlot <- "Random locus quantiles"

# Define quantile colours
quantileColours <- c("red", "navy")
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
# Subset chrs to only those within a given subgenome
genomeLetter <- unlist(strsplit(gsub("NLR_genes_in_", "",
                                     gsub("genome_genomewide", "", featureName)),
                                split = "_"))
if(length(genomeLetter) == 1) {
  chrs <- chrs[grepl(genomeLetter, chrs)]
}

# Load table of features grouped into quantiles
# by decreasing log2(libName/control)
mclapply(seq_along(pop_name), function(x) {
if(libName %in% c("cMMb", "cluster_members", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  featuresDF <- read.table(paste0(outDir[x], "features_", quantiles, "quantiles",
                                  "_by_", libName, "_of_",
                                  substring(featureName[1][1], first = 1, last = 9), "_in_",
                                  paste0(substring(featureName, first = 14, last = 20),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 22), "_", pop_name[x], ".txt"),
                           header = T, sep = "\t")
} else {
  featuresDF <- read.table(paste0(outDir[x], "features_", quantiles, "quantiles",
                                  "_by_", sub("_\\w+", "", libName), "_in_",
                                  region, "_of_",
                                  substring(featureName[1][1], first = 1, last = 9), "_in_",
                                  paste0(substring(featureName, first = 14, last = 20),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 22), "_", pop_name[x], ".txt"),
                           header = T, sep = "\t")
}

featuresDF <- featuresDF[featuresDF$NLR_quantile != "",]
rownames(featuresDF) <- as.character(1:dim(featuresDF)[1])

# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$NLR_quantile == paste0("Quantile ", k))
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
#set.seed(93750174)
set.seed(453838430)

# Divide features into random sets of equal number,
# with the same number of genes per chromosome as
# above-defined libName-defined feature quantiles
randomPCIndices <- lapply(1:quantiles, function(k) {
  randomPCIndicesk <- NULL
  for(i in 1:length(chrs)) {
    randomPCfeatureskChr <- selectRandomFeatures(features = featuresDF[featuresDF$seqnames == chrs[i],],
                                                 n = dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k) &
                                                                    featuresDF$seqnames == chrs[i],])[1])
    if(dim(randomPCfeatureskChr)[1] > 0) {
      randomPCIndicesk <- c(randomPCIndicesk, as.integer(rownames(randomPCfeatureskChr)))
    }
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
## Alternatively, randomly select without considering per-chromosome features
#randomPCIndices <- lapply(1:quantiles, function(k) {
#  randomPCfeatureskChr <- selectRandomFeatures(features = featuresDF,
#                                               n = dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),])[1])
#  as.integer(rownames(randomPCfeatureskChr))
#})

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
  mean(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),][,colnames(featuresDF) == orderingFactor,], na.rm = T)
})
featuresDF_quantileSD <- sapply(1:quantiles, function(k) {
  sd(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),][,colnames(featuresDF) == orderingFactor,], na.rm = T)
})
featuresDF_quantileSEM <- sapply(1:quantiles, function(k) {
  featuresDF_quantileSD[k] / sqrt( (dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),])[1] - 1) )
})
featuresDF_quantileCIlower <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] -
    ( qt(0.975, df = dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_quantileCIupper <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] +
    ( qt(0.975, df = dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_summary_stats <- data.frame(NLR_quantile = paste0("Quantile ", 1:quantiles),
                                       Mean = featuresDF_quantileMean,
                                       SD = featuresDF_quantileSD,
                                       SEM = featuresDF_quantileSEM,
                                       CIlower = featuresDF_quantileCIlower,
                                       CIupper = featuresDF_quantileCIupper,
                                       stringsAsFactors = F)

ranFeatsDF_randomMean <- sapply(1:quantiles, function(k) {
  mean(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),][,colnames(ranFeatsDF) == orderingFactor,], na.rm = T)
})
ranFeatsDF_randomSD <- sapply(1:quantiles, function(k) {
  sd(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),][,colnames(ranFeatsDF) == orderingFactor,], na.rm = T)
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
ranFeatsDF_summary_stats <- data.frame(random = paste0("Random ", 1:quantiles),
                                       Mean = ranFeatsDF_randomMean,
                                       SD = ranFeatsDF_randomSD,
                                       SEM = ranFeatsDF_randomSEM,
                                       CIlower = ranFeatsDF_randomCIlower,
                                       CIupper = ranFeatsDF_randomCIupper,
                                       stringsAsFactors = F)
summary_stats_min <- min(c(featuresDF_summary_stats$CIlower, ranFeatsDF_summary_stats$CIlower), na.rm = T)
summary_stats_max <- max(c(featuresDF_summary_stats$CIupper, ranFeatsDF_summary_stats$CIupper), na.rm = T)

#k <- 1
#dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k) &
#               !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor)]) &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] != 0 &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] < 0.001,])
#
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] < 0.0002,])
#dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k) &
##               !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor)]) &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] >= 0 &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] < 0.0002,])
#dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k) &
##               !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor)]) &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] >= 0.0002 &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] < 0.0002,])
#max(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k) &
##               !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor)]) &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] >= 0 &
#               featuresDF[,which(colnames(featuresDF) == orderingFactor)] < 0.0002,][,which(colnames(featuresDF) == orderingFactor)],
#    na.rm = T)
#
#dim(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k) &
#               ranFeatsDF[,which(colnames(ranFeatsDF) == orderingFactor)] >= 0 &
#               ranFeatsDF[,which(colnames(ranFeatsDF) == orderingFactor)] < 0.02,])

# Take top densityProp*100% of data to aid visualisation in density plots
featuresDF <- featuresDF[which(featuresDF[,which(colnames(featuresDF) == orderingFactor)] <=
                               quantile(featuresDF[,which(colnames(featuresDF) == orderingFactor)],
                                        probs = densityProp, na.rm = T)),]
#                               featuresDF[,which(colnames(featuresDF) == orderingFactor)] != 0,]
ranFeatsDF <- ranFeatsDF[which(ranFeatsDF[,which(colnames(ranFeatsDF) == orderingFactor)] <=
                               quantile(ranFeatsDF[,which(colnames(ranFeatsDF) == orderingFactor)],
                                        probs = densityProp, na.rm = T)),]
#                               ranFeatsDF[,which(colnames(ranFeatsDF) == orderingFactor)] != 0,]
xmin <- min(c(featuresDF[,which(colnames(featuresDF) == orderingFactor)]),
              na.rm = T)
xmax <- max(c(featuresDF[,which(colnames(featuresDF) == orderingFactor)]),
              na.rm = T)
minDensity <- 0
maxDensity <- max(density(featuresDF[featuresDF$NLR_quantile == "Quantile 2",][,which(colnames(featuresDF) == orderingFactor)],
                          na.rm = T)$y)+maxDensityPlus
maxDensity <- max(
  c(
    sapply(1:quantiles, function(k) {
#      max(c(max(density(ranLocsDF[ranLocsDF$random == paste0("Random ", k),][,which(colnames(featuresDF) == orderingFactor)],
#                        na.rm = T)$y),
      max(c(max(density(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),][,which(colnames(featuresDF) == orderingFactor)],
                        na.rm = T)$y),
            max(density(ranFeatsDF[ranFeatsDF$random == paste0("Random ", k),][,which(colnames(featuresDF) == orderingFactor)],
                        na.rm = T)$y)))
     })
   )
)+maxDensityPlus

# Define legend labels
legendLabs_feature <- lapply(1:quantiles, function(x) {
  grobTree(textGrob(bquote(.(paste0("Quantile ", 1:quantiles)[x])),
                    x = legendLabX, y = 0.90-((x-1)*0.07), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 22)))
})
legendLabs_ranFeat <- lapply(1:quantiles, function(x) {
  grobTree(textGrob(bquote(.(paste0("Random ", 1:quantiles)[x])),
                    x = legendLabX, y = 0.90-((x-1)*0.07), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 22)))
})

# Population genetics statistic density plot function
popgen_stats_plotFun <- function(lociDF,
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
                     labels = function(x) sprintf(xDec, x)) +
  scale_y_continuous(limits = c(minDensity, maxDensity),
                     labels = function(x) sprintf(yDec, x)) +
  labs(x = parameterLab,
       y = "Density") +
  annotation_custom(legendLabs[[1]]) +
  annotation_custom(legendLabs[[2]]) +
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
popgen_stats_meanCIs <- function(dataFrame,
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
                     labels = function(x) sprintf(yDec2, x)) +
#  scale_x_discrete(breaks = as.vector(dataFrame$NLR_quantile),
#                   labels = as.vector(dataFrame$NLR_quantile)) +
  labs(x = "",
       y = parameterLab) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.text.x = element_markdown(size = 22, colour = quantileColours, hjust = 1.0, vjust = 1.0, angle = 45),
        axis.title = element_text(size = 26, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.1,0.3),"cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot)))
}

ggObjGA_feature <- popgen_stats_plotFun(lociDF = featuresDF[grepl("Quantile ", featuresDF$NLR_quantile),],
                                        parameter = orderingFactor,
                                        parameterLab = bquote(.(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
                                        featureGroup = "NLR_quantile", 
                                        featureNamePlot = featureNamePlot,
                                        legendLabs = legendLabs_feature,
                                        quantileColours = quantileColours
                                       )
ggObjGA_ranFeat <- popgen_stats_plotFun(lociDF = ranFeatsDF[grepl("Random ", ranFeatsDF$random),],
                                        parameter = orderingFactor,
                                        parameterLab = bquote(.(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
                                        featureGroup = "random", 
                                        featureNamePlot = ranFeatNamePlot,
                                        legendLabs = legendLabs_ranFeat,
                                        quantileColours = quantileColours
                                       )
ggObjGA_feature_mean <- popgen_stats_meanCIs(dataFrame = featuresDF_summary_stats,
                                             parameterLab = bquote(.(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
                                             featureGroup = "NLR_quantile",
                                             featureNamePlot = featureNamePlot,
                                             quantileColours = quantileColours
                                            )
ggObjGA_ranFeat_mean <- popgen_stats_meanCIs(dataFrame = ranFeatsDF_summary_stats,
                                             parameterLab = bquote(.(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
                                             featureGroup = "random",
                                             featureNamePlot = ranFeatNamePlot,
                                             quantileColours = quantileColours
                                            )
ggObjGA_combined <- grid.arrange(ggObjGA_feature,
                                 ggObjGA_feature_mean,
                                 ggObjGA_ranFeat,
                                 ggObjGA_ranFeat_mean,
                                 ncol = 2, as.table = F)
if(libName %in% c("cMMb", "cluster_members", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir[x],
                orderingFactor, "_densityProp", densityProp, "_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 9), "_in_",
                paste0(substring(featureName, first = 14, last = 20),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 22), "_", pop_name[x], "_v250620.pdf"),
         plot = ggObjGA_combined,
         height = 13, width = 14)
} else {  
  ggsave(paste0(plotDir[x],
                orderingFactor, "_densityProp", densityProp, "_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 9), "_in_",
                paste0(substring(featureName, first = 14, last = 20),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 22), "_", pop_name[x], "_v250620.pdf"),
         plot = ggObjGA_combined,
         height = 13, width = 14)
}
}, mc.cores = length(pop_name), mc.preschedule = F)
