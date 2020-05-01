#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage profiles with 95% CIs around peak quantiles

# Usage:
# csmit -m 100G -c 15 "/applications/R/R-3.5.0/bin/Rscript quantile_peaks_avgProfileRibbon_cMMb_TEsuperfams.R DMC1_Rep1_ChIP DMC1 'Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin' cMMb 4 both 400 2000 2kb 20 '0.02,0.96'"

#libName <- "DMC1_Rep1_ChIP"
#dirName <- "DMC1"
#featureName <- unlist(strsplit("Agenome_euchromatin,Bgenome_euchromatin,Dgenome_euchromatin",
#                               split = ","))
#orderingFactor <- "cMMb"
#quantiles <- 4
#align <- "both"
#bodyLength <- 400
#upstream <- 2000
#downstream <- 2000
#flankName <- "2 kb"
#binSize <- 20
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
#                                        split = ",")))
args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
orderingFactor <- args[4]
quantiles <- as.numeric(args[5])
align <- args[6]
bodyLength <- as.numeric(args[7])
upstream <- as.numeric(args[8])
downstream <- as.numeric(args[8])
flankName <- args[9]
binSize <- as.numeric(args[10])
legendPos <- as.numeric(unlist(strsplit(args[11],
                                        split = ",")))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

outDir <- paste0("quantiles_by_", orderingFactor, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Define plot titles
if(orderingFactor == "cMMb") {
  featureNamePlot <- paste0(substr(orderingFactor, start = 1, stop = 2), "/",
                            substr(orderingFactor, start = 3, stop = 4), " ",
                            sub("_\\w+", "", libName), " peaks")
}
ranFeatNamePlot <- paste0("Random ",
                          sub("_\\w+", "", libName), " peaks")
ranLocNamePlot <- "Random loci"

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]

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

# Divide features into random sets of equal number,
# with the same number of peaks per chromosome as
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
                  dim(featuresDF[quantileIndices[[k]],][featuresDF[quantileIndices[[k]],]$seqnames == chrs[x],])))    {
      stop("Quantile features and random features do not consist of the same number of features per chromosome")
    }
  })
})


# TE superfams
superfamCodes <- c("RLG",
                   "RLC",
                   "RLX",
                   "RIX",
                   "SIX",
                   "DTC",
                   "DTM",
                   "DTX",
                   "DTH",
                   "DMI",
                   "DTT",
                   "DXX",
                   "DTA",
                   "DHH",
                   "XXX")
superfamNames <- c("Gypsy_LTR",
                   "Copia_LTR",
                   "Unclassified_LTR",
                   "LINE",
                   "SINE",
                   "CACTA",
                   "Mutator",
                   "Unclassified_with_TIRs",
                   "Harbinger",
                   "MITE",
                   "Mariner",
                   "Unclassified_class_2",
                   "hAT",
                   "Helitrons",
                   "Unclassified_repeats")
superfamNamesPlot <- c("Gypsy LTR",
                       "Copia LTR",
                       "Unclassified LTR",
                       "LINE",
                       "SINE",
                       "CACTA",
                       "Mutator",
                       "Unclassified with TIRs",
                       "Harbinger",
                       "MITE",
                       "Mariner",
                       "Unclassified class 2",
                       "hAT",
                       "Helitrons",
                       "Unclassified repeats")

# feature
superfam_featureMats <- mclapply(seq_along(superfamNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/DMC1peakProfiles/matrices/",
                                superfamNames[x], "_", superfamCodes[x],
                                "_around_", dirName, "_peaks_in_", featureName[y],
                                "_matrix_bin", binSize, "bp_flank", sub(" ", "", flankName), ".tab"),
                         header = T))
  })
}, mc.cores = length(superfamNames))
# If features from multiple subgenomes and/or compartments are to be analysed,
# concatenate the corresponding feature coverage matrices
superfam_featureMats <- mclapply(seq_along(superfam_featureMats), function(x) {
  if(length(featureName) > 1) {
    do.call(rbind, superfam_featureMats[[x]])
  } else {
    superfam_featureMats[[x]][[1]]
  }
}, mc.cores = length(superfam_featureMats))

# ranLoc
superfam_ranLocMats <- mclapply(seq_along(superfamNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/DMC1peakProfiles/matrices/",
                                superfamNames[x], "_", superfamCodes[x],
                                "_around_", dirName, "_peaks_in_", featureName[y],
                                "_ranLoc_matrix_bin", binSize, "bp_flank", sub(" ", "", flankName), ".tab"),
                         header = T))
  })
}, mc.cores = length(superfamNames))
# If features from multiple subgenomes and/or compartments are to be analysed,
# concatenate the corresponding feature coverage matrices
superfam_ranLocMats <- mclapply(seq_along(superfam_ranLocMats), function(x) {
  if(length(featureName) > 1) {
    do.call(rbind, superfam_ranLocMats[[x]])
  } else {
    superfam_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(superfam_ranLocMats))

# Add column names
for(x in seq_along(superfam_featureMats)) {
  colnames(superfam_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(superfam_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
superfam_mats_quantiles <- mclapply(seq_along(superfam_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         superfam_ranLocMats[[x]][quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(superfam_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_superfam <- mclapply(seq_along(superfam_mats_quantiles), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(superfam_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(superfam_mats_quantiles[[x]][[y]][[k]]),
                 t(superfam_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(superfam_mats_quantiles)/3)

# Convert into tidy data.frame (long format)
tidyDFfeature_list_superfam  <- mclapply(seq_along(wideDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(superfam_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_superfam[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_superfam)/3)

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_superfam)) {
  for(y in seq_along(superfam_mats_quantiles[[x]])) {
    for(k in seq_along(superfam_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_superfam[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_superfam[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_superfam  <- mclapply(seq_along(tidyDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(superfam_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_superfam[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_superfam)/3)

for(x in seq_along(summaryDFfeature_list_superfam)) {
  for(y in seq_along(superfam_mats_quantiles[[x]])) {
    for(k in seq_along(superfam_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_superfam[[x]][[y]][[k]]$window))
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_superfam[[x]][[y]][[k]])[1])
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$sem <- summaryDFfeature_list_superfam[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_superfam[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_superfam[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_superfam[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_superfam[[x]][[y]][[k]]$sem
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_superfam[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_superfam[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_superfam[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_superfam)) {
  # feature quantiles
  names(summaryDFfeature_list_superfam[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_superfam[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_superfam[[x]][[3]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_superfam into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_superfam  <- mclapply(seq_along(summaryDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_superfam[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_superfam))
for(x in seq_along(summaryDFfeature_superfam)) {
  # feature quantiles
  summaryDFfeature_superfam[[x]][[1]]$quantile <- factor(summaryDFfeature_superfam[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_superfam[[x]][[2]]$quantile <- factor(summaryDFfeature_superfam[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_superfam[[x]][[3]]$quantile <- factor(summaryDFfeature_superfam[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[3]]))
}

# Define y-axis limits
ymin_list_superfam <- lapply(seq_along(summaryDFfeature_superfam), function(x) {
  min(c(summaryDFfeature_superfam[[x]][[1]]$CI_lower,
        summaryDFfeature_superfam[[x]][[2]]$CI_lower,
        summaryDFfeature_superfam[[x]][[3]]$CI_lower))
})
ymax_list_superfam <- lapply(seq_along(summaryDFfeature_superfam), function(x) {
  max(c(summaryDFfeature_superfam[[x]][[1]]$CI_upper,
        summaryDFfeature_superfam[[x]][[2]]$CI_upper,
        summaryDFfeature_superfam[[x]][[3]]$CI_upper))
})

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_feature[[1]]) +
  annotation_custom(legendLabs_feature[[2]]) +
  annotation_custom(legendLabs_feature[[3]]) +
  annotation_custom(legendLabs_feature[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranFeat
ggObj2_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranLoc
ggObj3_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              "Start",
                              "End",
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_superfam,
                                           ggObj2_combined_superfam,
                                           ggObj3_combined_superfam
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(superfamNamesPlot)),
                                                       (length(c(superfamNamesPlot))+1):(length(c(superfamNamesPlot))*2),
                                                       ((length(c(superfamNamesPlot))*2)+1):(length(c(superfamNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "TEsuperfam_avgProfiles_around_", quantiles, "quantiles",
               "_by_", orderingFactor,
               "_of_", libName, "_peaks_in_",
               paste0(featureName,
                      collapse = "_"), "_v300420.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(superfamNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   superfam_featureMats, superfam_ranLocMats,
   superfam_mats_quantiles,
   wideDFfeature_list_superfam,
   tidyDFfeature_list_superfam,
   summaryDFfeature_list_superfam,
   summaryDFfeature_superfam
  ) 
gc()
#####
