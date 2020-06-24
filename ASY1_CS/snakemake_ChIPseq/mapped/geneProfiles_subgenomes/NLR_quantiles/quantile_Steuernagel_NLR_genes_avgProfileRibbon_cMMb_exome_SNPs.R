#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage profiles with 95% CIs around
# gene quantiles

# Usage:
# /applications/R/R-3.5.0/bin/Rscript quantile_Steuernagel_NLR_genes_avgProfileRibbon_cMMb_exome_SNPs.R cMMb both 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' 3500 2000 2kb '2 kb' 20 20bp genes 4 100kb 1 '0.02,0.96'

#libName <- "cMMb"
#align <- "both"
#featureName <- unlist(strsplit("NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide",
#                               split = ","))
#bodyLength <- 3500
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "genes"
#quantiles <- 4
#winName <- "100kb"
#minMarkerDist <- 1
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
align <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
bodyLength <- as.numeric(args[4])
upstream <- as.numeric(args[5])
downstream <- as.numeric(args[5])
flankName <- args[6]
flankNamePlot <- args[7]
binSize <- as.numeric(args[8])
binName <- args[9]
region <- args[10]
quantiles <- as.numeric(args[11])
winName <- args[12]
minMarkerDist <-  as.numeric(args[13])
legendPos <- as.numeric(unlist(strsplit(args[14],
                                        split = ",")))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

outDir <- paste0("quantiles_by_", libName, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Define plot titles
if(libName == "cMMb") {
  featureNamePlot <- paste0("cM/Mb ",
                            substr(featureName[1], start = 1, stop = 3),
                            " quantiles")
} else if(libName == "cluster_members") {
  featureNamePlot <- paste0(substr(featureName[1], start = 1, stop = 3),
                            " cluster size",
                            " quantiles")
} else {
  featureNamePlot <- paste0(libName, " ",
                            substr(featureName[1], start = 1, stop = 3),
                            " quantiles (", region, ")")
}
ranFeatNamePlot <- paste0("Random ",
                          substr(featureName[1], start = 1, stop = 3),
                          " quantiles")
#ranLocNamePlot <- "Random locus quantiles"

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")

# Define feature start and end labels for plotting
if(grepl("genes", featureName)) {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]

# Load table of features grouped into quantiles
featuresDF <- read.table(paste0(outDir,
                                "/WesternEurope/features_", quantiles, "quantiles_by_",
                                libName, "_of_",
                                substring(featureName[1][1], first = 1, last = 9), "_in_",
                                paste0(substring(featureName, first = 14, last = 20),
                                       collapse = "_"), "_",
                                substring(featureName[1][1], first = 22), "_WesternEurope.txt"),
                         header = T, sep = "\t")

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
features <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                    "IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                    paste0(substring(featureName[x], first = 14, last = 20),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 22), ".gff3"),
             header = F)
})
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature data.frames
if(length(featureName) == 3) {
 features <- do.call(rbind, features)
} else {
 features <- features[[1]]
}

# Remove NLRs that have not been assigned to a quantile
# due to missing libName data
featuresDF <- featuresDF[featuresDF$NLR_quantile != "",]
rownames(featuresDF) <- as.character(1:dim(featuresDF)[1])

# Get row indices for each feature quantile
# The coverage matrix contains coverage profiles for all genes,
# not just NLR-encoding genes, so it's necessary to get row indices
# corresponding to NLRs in each quantile
quantileNLRIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$NLR_quantile == paste0("Quantile ", k))
})
quantileIndices <- lapply(1:quantiles, function(k) {
  which(as.character(features$V9) %in% as.character(featuresDF[quantileNLRIndices[[k]],]$featureID))
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
randomPCNLRIndices <- lapply(1:quantiles, function(k) {
  randomPCIndicesk <- NULL
  for(x in 1:length(chrs)) {
    randomPCfeatureskChr <- selectRandomFeatures(features = featuresDF[featuresDF$seqnames == chrs[x],],
                                                 n = dim(featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k) &
                                                                    featuresDF$seqnames == chrs[x],])[1])
    if(dim(randomPCfeatureskChr)[1] > 0) {
      randomPCIndicesk <- c(randomPCIndicesk, as.integer(rownames(randomPCfeatureskChr)))
    }
  }
  randomPCIndicesk
})
# Confirm per-chromosome feature numbers are the same for quantiles and random groupings
lapply(seq_along(1:quantiles), function(k) {
  print(k)
  sapply(seq_along(chrs), function(x) {
    print(x)
    if(!identical(dim(featuresDF[randomPCNLRIndices[[k]],][featuresDF[randomPCNLRIndices[[k]],]$seqnames == chrs[x],]),
                  dim(featuresDF[quantileNLRIndices[[k]],][featuresDF[quantileNLRIndices[[k]],]$seqnames == chrs[x],])))     {
      stop("Quantile features and random features do not consist of the same number of features per chromosome")
    }
  })
})
# The coverage matrix contains coverage profiles for all genes,
# not just NLR-encoding genes, so it's necessary to get row indices
# corresponding to NLRs in each quantile
randomPCIndices <- lapply(1:quantiles, function(k) {
  which(as.character(features$V9) %in% as.character(featuresDF[randomPCNLRIndices[[k]],]$featureID))
})

# exome SNPclasses
SNPclassNames <- c(
                   "all",
                   "transition",
                   "transversion"
                  )
SNPclassNamesPlot <- c(
                       "Exome SNPs",
                       "Transitions",
                       "Transversions"
                      )

# feature
SNPclass_featureMats <- mclapply(seq_along(SNPclassNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/exome_",
                                SNPclassNames[x],
                                "_SNPs_around_", featureName[y],
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(SNPclassNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
SNPclass_featureMats <- mclapply(seq_along(SNPclass_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, SNPclass_featureMats[[x]])
  } else {
    SNPclass_featureMats[[x]][[1]]
  }
}, mc.cores = length(SNPclass_featureMats))

# ranLoc
SNPclass_ranLocMats <- mclapply(seq_along(SNPclassNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/exome_",
                                SNPclassNames[x],
                                "_SNPs_around_", featureName[y],
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(SNPclassNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
SNPclass_ranLocMats <- mclapply(seq_along(SNPclass_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, SNPclass_ranLocMats[[x]])
  } else {
    SNPclass_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(SNPclass_ranLocMats))

# Add column names
for(x in seq_along(SNPclass_featureMats)) {
  colnames(SNPclass_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(SNPclass_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
SNPclass_mats_quantiles <- mclapply(seq_along(SNPclass_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         SNPclass_ranLocMats[[x]][quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(SNPclass_featureMats))


# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_SNPclass <- mclapply(seq_along(SNPclass_mats_quantiles), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(SNPclass_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(SNPclass_mats_quantiles[[x]][[y]][[k]]),
                 t(SNPclass_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(SNPclass_mats_quantiles)/2)

# Convert into tidy data.frame (long format)
tidyDFfeature_list_SNPclass  <- mclapply(seq_along(wideDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(SNPclass_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_SNPclass[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_SNPclass)/2)

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_SNPclass)) {
  for(y in seq_along(SNPclass_mats_quantiles[[x]])) {
    for(k in seq_along(SNPclass_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_SNPclass[[x]][[y]][[k]]$window))
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
summaryDFfeature_list_SNPclass  <- mclapply(seq_along(tidyDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(SNPclass_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_SNPclass[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_SNPclass)/2)

for(x in seq_along(summaryDFfeature_list_SNPclass)) {
  for(y in seq_along(SNPclass_mats_quantiles[[x]])) {
    for(k in seq_along(SNPclass_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_SNPclass[[x]][[y]][[k]]$window))
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_SNPclass[[x]][[y]][[k]])[1])
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sem <- summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sem
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_SNPclass)) {
  # feature quantiles
  names(summaryDFfeature_list_SNPclass[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_SNPclass[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_SNPclass[[x]][[3]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_SNPclass into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_SNPclass  <- mclapply(seq_along(summaryDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_SNPclass[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_SNPclass))
for(x in seq_along(summaryDFfeature_SNPclass)) {
  # feature quantiles
  summaryDFfeature_SNPclass[[x]][[1]]$NLR_quantile <- factor(summaryDFfeature_SNPclass[[x]][[1]]$NLR_quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_SNPclass[[x]][[2]]$NLR_quantile <- factor(summaryDFfeature_SNPclass[[x]][[2]]$NLR_quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_SNPclass[[x]][[3]]$NLR_quantile <- factor(summaryDFfeature_SNPclass[[x]][[3]]$NLR_quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[3]]))
}

# Define y-axis limits
ymin_list_SNPclass <- lapply(seq_along(summaryDFfeature_SNPclass), function(x) {
  min(c(summaryDFfeature_SNPclass[[x]][[1]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[2]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[3]]$CI_lower))
})
ymax_list_SNPclass <- lapply(seq_along(summaryDFfeature_SNPclass), function(x) {
  max(c(summaryDFfeature_SNPclass[[x]][[1]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[2]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[3]]$CI_upper))
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
ggObj1_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[1]]
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
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
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
}, mc.cores = length(SNPclassNamesPlot))

## ranFeat
ggObj2_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[2]]
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
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
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
}, mc.cores = length(SNPclassNamesPlot))

## ranLoc
ggObj3_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[3]]
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
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
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
}, mc.cores = length(SNPclassNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_SNPclass,
                                           ggObj2_combined_SNPclass,
                                           ggObj3_combined_SNPclass
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(SNPclassNamesPlot)),
                                                       (length(c(SNPclassNamesPlot))+1):(length(c(SNPclassNamesPlot))*2),
                                                       ((length(c(SNPclassNamesPlot))*2)+1):(length(c(SNPclassNamesPlot))*3)
                                                      ))
ggsave(paste0(plotDir,
              "exomeSNPclass_avgProfiles_around_", quantiles, "quantiles",
              "_by_", libName, "_of_",
              substring(featureName[1][1], first = 1, last = 5), "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18), "_v090620.pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(SNPclassNamesPlot)), width = 21, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   SNPclass_featureMats, SNPclass_ranLocMats,
   SNPclass_mats_quantiles,
   wideDFfeature_list_SNPclass,
   tidyDFfeature_list_SNPclass,
   summaryDFfeature_list_SNPclass,
   summaryDFfeature_SNPclass
  ) 
gc()
