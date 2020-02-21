#!/applications/R/R-3.5.0/bin/Rscript

# Compare population genetics statistics (means) for GO-term-annotated genes in a given ASY1 quantile
# with those for randomSets sets of randomly selected genes not in that ASY1 quantile, and
# with those for randomSets sets of randomly selected GO-term-annotated genes not in that ASY1 quantile

# Usage:
# /applications/R/R-3.5.0/bin/Rscript quantile_GOann_genes_popgen_stats_permTest.R ASY1_CS_Rep1_ChIP ASY1_CS 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' 'Defense_response_genes' '0006952' promoters 1 4 TajimaD 'Tajima D' 10000 0.0001

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#featureNamePlot <- "Defense_response_genes"
#GO_ID <- "0006952"
#region <- "promoters"
#quantileNo <- 1
#quantiles <- 4
#orderingFactor <- "TajimaD"
#orderingFactorName <- bquote("Tajima's" ~ italic("D"))
#orderingFactor <- "RozasR2"
#orderingFactorName <- bquote("Rozas'" ~ italic("R")[2])
#randomSets <- 10000
#minPval <- 0.0001

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
featureNamePlot <- args[4]
GO_ID <- args[5]
region <- args[6]
quantiles <- as.numeric(args[7])
orderingFactor <- args[8]
orderingFactorName <- args[9]
randomSets <- as.numeric(args[10])
minPval <- as.numeric(args[11])

library(parallel)
library(plotrix)
#library(tidyr)
#library(dplyr)
#library(ggplot2)
#library(ggbeeswarm)
#library(ggthemes)
#library(grid)
#library(gridExtra)
#library(extrafont)

pop_name <- c("NorthAfrica",
              "SubSaharanAfrica",
              "WesternEurope",
              "EasternEurope",
              "MiddleEast",
              "FormerSU",
              "CentralAsia",
              "SouthAsia",
              "EastAsia",
              "NorthAmerica",
              "CentralAmerica",
              "SouthAmerica",
              "Oceania")
pop_name_plot <- c("North Africa",
                   "Sub-Saharan Africa",
                   "Western Europe",
                   "Eastern Europe",
                   "Middle East",
                   "Former SU",
                   "Central Asia",
                   "South Asia",
                   "East Asia",
                   "North America",
                   "Central America",
                   "South America",
                   "Oceania")

outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                 "_in_", region, "/")
outDir <- sapply(seq_along(pop_name), function(x) {
  paste0(outDir, pop_name[x], "/")
})
plotDir <- paste0(outDir, "plots/")


# Load IDs of genes annotated with GO_ID in quantile quantileNo
IDs <- as.character(read.table(paste0("quantiles_by_log2_", libName, "_control_in_", region,
                                      "/GO/featureIDs_quantile", quantileNo, "_of_", quantiles,
                                      "_by_log2_", libName, "_control_in_", region, "_of_",
                                      substring(featureName[1][1], first = 1, last = 5), "_in_",
                                      paste0(substring(featureName, first = 10, last = 16),
                                             collapse = "_"), "_",
                                      substring(featureName[1][1], first = 18),
                                      "_GO_BP/featureIDs_quantile", quantileNo, "_of_", quantiles,
                                      "_by_log2_", libName, "_control_in_", region, "_of_",
                                      substring(featureName[1][1], first = 1, last = 5), "_in_",
                                      paste0(substring(featureName, first = 10, last = 16),
                                             collapse = "_"), "_",
                                      substring(featureName[1][1], first = 18),
                                      "_GO_BP_enrichment_GO:", GO_ID, ".txt"),
                               colClasses = c("NULL", NA), header = F)$V2)
IDs <- unlist(strsplit(x = IDs,
                       split = ","))

# Load table of functional annotations
# (note this is old v1.0 and not new v1.1 gene annotation)
anno <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned.TAB",
                   header = T, quote = "\"", sep = "\t", check.names = F)

# Replace gene model ID decimal suffix (e.g., ".1")
anno$`Gene-ID` <- sub(pattern = "\\.\\d+", replacement = "",
                      x = anno$`Gene-ID`)
# Replace "1G" with "2G" in gene IDs for consistency with v1.1
anno$`Gene-ID` <- sub(pattern = "1G", replacement = "2G",
                      x = anno$`Gene-ID`)
# Get subset corresponding to genes annotated with enriched GO_ID, and not in quantile quantileNo
annoGOIDs <- anno[which(grepl(pattern = GO_ID, x = anno$`GO-IDs-via-Interpro`,
                              fixed = T)),]
annoGOIDs <- annoGOIDs[!(annoGOIDs$`Gene-ID` %in% IDs),]$`Gene-ID`

# Load table of features grouped into quantiles
# by decreasing log2(libName/control) in region
#mclapply(seq_along(pop_name), function(x) {
for(x in seq_along(pop_name)) {
  featuresDF <- read.table(paste0(outDir[x], "features_", quantiles, "quantiles",
                                  "_by_", sub("_\\w+", "", libName), "_in_",
                                  region, "_of_",
                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
                                  paste0(substring(featureName, first = 10, last = 16),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                           header = T, sep = "\t", stringsAsFactors = F)
  featuresDF$featureID <- sub("\\.\\d+", "", featuresDF$featureID)
  # Retain rows that contain non-NA values of orderingFactor
  featuresDF <- featuresDF[!is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor)]),]
  
  # Get orderingFactor values for IDs
  featuresDF_IDsDF <- featuresDF[featuresDF$featureID %in% IDs,]
  #featuresDF_nonIDsDF <- featuresDF[!(featuresDF$featureID %in% IDs),]
  featuresDF_nonIDsDF <- featuresDF[featuresDF$quantile != paste0("Quantile ", quantileNo),]
  featuresDF_annoGOIDsDF <- featuresDF[featuresDF$featureID %in% annoGOIDs,]
  featuresDF_IDs <- featuresDF_IDsDF$featureID
  featuresDF_nonIDs <- featuresDF_nonIDsDF$featureID
  featuresDF_annoGOIDs <- featuresDF_annoGOIDsDF$featureID
  
  # Calculate mean orderingFactor value for genes annotated with GO_ID in quantile quantileNo
  IDs_mean <- mean(featuresDF_IDsDF[,which(colnames(featuresDF_IDsDF) ==
                                           orderingFactor)],
                   na.rm = T)
  
  # Function to randomly select n feature IDs not present in IDs (from featuresDF_nonIDs or featuresDF_annoGOIDs)
  ran_nonIDs_select <- function(nonIDsChr, n) {
    sample(x = nonIDsChr,
           size = n,
           replace = FALSE)
  }
  
  # Define seed so that random selections are reproducible
  set.seed(453838430)
  
  # Apply ran_nonIDs_select() function on a per-chromosome basis to generate randomSets random sets
  # of IDs from featuresDF_nonIDs for which to calculate randomSets mean orderingFactor values
  ran_nonIDs_perm_means <- unlist(mclapply(1:randomSets, function(y) {
    ran_nonIDs <- NULL
    for(i in 1:length(unique(featuresDF_IDsDF$seqnames))) {
      IDsChr <- featuresDF_IDs[grepl(paste0("TraesCS",
                                            sub("chr", "", unique(featuresDF_IDsDF$seqnames))[i]),
                                     featuresDF_IDs)]
      nonIDsChr <- featuresDF_nonIDs[grepl(paste0("TraesCS",
                                                  sub("chr", "", unique(featuresDF_IDsDF$seqnames))[i]),
                                           featuresDF_nonIDs)]
      ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                         n = length(IDsChr))
      ran_nonIDs <- c(ran_nonIDs, ran_nonIDsChr)
    }
    return(mean(featuresDF_nonIDsDF[featuresDF_nonIDsDF$featureID %in%
                                    ran_nonIDs,][,which(colnames(featuresDF_nonIDsDF) ==
                                                        orderingFactor)],
                na.rm = T))
  }, mc.cores = detectCores(), mc.preschedule = T))
  
  # Define seed so that random selections are reproducible
  set.seed(453838430)
  
  # Apply ran_nonIDs_select() function to generate randomSets random sets
  # of IDs from featuresDF_annoGOIDs for which to calculate randomSets mean orderingFactor values
  ran_annoGOIDs_perm_means <- unlist(mclapply(1:randomSets, function(y) {
    ran_annoGOIDs <- ran_nonIDs_select(nonIDsChr = featuresDF_annoGOIDs,
                                       n = length(featuresDF_IDs))
    return(mean(featuresDF_annoGOIDsDF[featuresDF_annoGOIDsDF$featureID %in%
                                       ran_annoGOIDs,][,which(colnames(featuresDF_annoGOIDsDF) ==
                                                              orderingFactor)],
                na.rm = T))
  }, mc.cores = detectCores(), mc.preschedule = T))
  
  # Set class for permutation test results object
  setClass("permTest_popgen_stats",
           representation(alternative = "character",
                          alpha0.05 = "numeric",
                          pval = "numeric",
                          observed = "numeric",
                          permuted = "numeric",
                          expected = "numeric",
                          log2obsexp = "numeric",
                          log2alpha = "numeric"))
  
  # Disable scientific notation (e.g., 0.0001 rather than 1e-04)
  options(scipen = 100)
  
  # Determine whether mean orderingFactor values for genes annotated with GO_ID in quantile quantileNo
  # are lower than or higher than those for genes randomly selected from featuresDF_nonIDs and featuresDF_annoGOIDs
  IDs_lessThan_nonIDs_Bool <- IDs_mean < ran_nonIDs_perm_means
  IDs_moreThan_nonIDs_Bool <- IDs_mean > ran_nonIDs_perm_means
  
  # Calculate P-values and significance levels
  if(IDs_mean < mean(ran_nonIDs_perm_means)) {
    nonIDs_pval <- 1 - ( sum(IDs_lessThan_nonIDs_Bool) / length(ran_nonIDs_perm_means) )
    if(nonIDs_pval == 0) {
      nonIDs_pval <- minPval
    }
    nonIDs_MoreOrLessThanRandom <- "LessThanRandom"
    nonIDs_alpha0.05 <- quantile(ran_nonIDs_perm_means, probs = 0.05)[[1]]
  } else {
    nonIDs_pval <- 1 - ( sum(IDs_moreThan_nonIDs_Bool) / length(ran_nonIDs_perm_means) )
    if(nonIDs_pval == 0) {
      nonIDs_pval <- minPval
    }
    nonIDs_MoreOrLessThanRandom <- "MoreThanRandom"
    nonIDs_alpha0.05 <- quantile(ran_nonIDs_perm_means, probs = 0.95)[[1]]
  }
  
  # Create permutation test results object
  nonIDs_permTestResults <- new("permTest_popgen_stats",
                                alternative = nonIDs_MoreOrLessThanRandom,
                                alpha0.05 = nonIDs_alpha0.05,
                                pval = nonIDs_pval,
                                observed = IDs_mean,
                                permuted = ran_nonIDs_perm_means,
                                expected = mean(ran_nonIDs_perm_means),
                                log2obsexp = log2(IDs_mean / mean(ran_nonIDs_perm_means)),
                                log2alpha = log2(nonIDs_alpha0.05 / mean(ran_nonIDs_perm_means)))
  save(nonIDs_permTestResults,
       file = paste0(outDir[x],
                     orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", featureNamePlot,
                     "_in_quantile", quantileNo, "_of_", quantiles,
                     "_by_log2_", libName, "_control_in_", region, "_of_",
                     substring(featureName[1][1], first = 1, last = 5), "_in_",
                     paste0(substring(featureName, first = 10, last = 16),
                            collapse = "_"), "_",
                     substring(featureName[1][1], first = 18),
                     "_ann_with_GO_BP_enrichment_GO:", GO_ID, ".RData"))

  # Generate histogram
  pdf(paste0(plotDir[x],
             orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", featureNamePlot,
             "_in_quantile", quantileNo, "_of_", quantiles,
             "_by_log2_", libName, "_control_in_", region, "_of_",
             substring(featureName[1][1], first = 1, last = 5), "_in_",
             paste0(substring(featureName, first = 10, last = 16),
                    collapse = "_"), "_",
             substring(featureName[1][1], first = 18),
             "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_hist.pdf"), 
             height = 4.5, width = 5)
  par(mar = c(3.1, 3.1, 4.1, 1.1),
      mgp = c(1.85, 0.75, 0))
  ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
  #options(scipen = 100)
  # Calculate max density
  maxDensityPlus <- max(density(nonIDs_permTestResults@permuted)$y)*1.2
  if(orderingFactor %in% c("TajimaD", "FuLiF", "FuLiD")) {
    if(nonIDs_permTestResults@alternative == "MoreThanRandom") {
      xlim <- c(pmin(-1, min(nonIDs_permTestResults@permuted)*1.2),
                pmax(nonIDs_permTestResults@observed/1.2, nonIDs_permTestResults@alpha0.05/1.2))
      textX1 <- quantile(xlim, 0.25)[[1]]
#      textX1 <- min(nonIDs_permTestResults@permuted)*1.15
    } else {
      xlim <- c(pmin(-1, nonIDs_permTestResults@observed*1.2),
                max(nonIDs_permTestResults@permuted)/1.2)
      textX1 <- quantile(xlim, 0.75)[[1]]
#      textX1 <- min(nonIDs_permTestResults@permuted)*1.15
    }
  } else {
    if(nonIDs_permTestResults@alternative == "MoreThanRandom") {
      xlim <- c(pmin(-1, min(nonIDs_permTestResults@permuted)/1.2),
                pmax(nonIDs_permTestResults@observed*1.2, nonIDs_permTestResults@alpha0.05*1.2))
      textX1 <- quantile(xlim, 0.25)[[1]]
#      textX1 <- min(nonIDs_permTestResults@permuted)/1.15
    } else {
      xlim <- c(pmin(-1, nonIDs_permTestResults@observed/1.2),
                max(nonIDs_permTestResults@permuted)*1.2)
      textX1 <- quantile(xlim, 0.75)[[1]]
#      textX1 <- min(nonIDs_permTestResults@permuted)/1.15
    }
  }
  hist(nonIDs_permTestResults@permuted,
       freq = FALSE,
       col = "dodgerblue",
       border = NA,
       lwd = 2,
       xlim = c(pretty(xlim)[1],
                pretty(xlim)[length(pretty(xlim))]),
       ylim = c(0,
                maxDensityPlus),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "",
       axes = FALSE)
  axis(side = 2,
       at = pretty(density(nonIDs_permTestResults@permuted)$y),
       lwd = 2)
  mtext(side = 2,
        text = "Density",
        line = 1.85)
  axis(side = 1,
       at = pretty(xlim),
       lwd = 2)
  mtext(side = 1,
        text = bquote("Mean" ~ .(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
        line = 1.85)
  titleText <- list(bquote(.(gsub("_", " ", featureNamePlot)) ~
                           "in" ~ .(sub("_\\w+", "", libName)) ~ "Quantile" ~ .(as.character(quantileNo)) ~
                           "(" * .(region) * ") in" ~
                           .(paste0(substring(featureName, first = 10, last = 16),
                                    collapse = " ")) ~ 
                           .(substring(featureName[1][1], first = 18))),
                    bquote(italic("P")*" = "*
                           .(as.character(round(nonIDs_permTestResults@pval,
                                                digits = 6)))),
                    bquote("Permutations = "*.(prettyNum(length(nonIDs_permTestResults@permuted),
                                                         big.mark = ",",
                                                         trim = T))))
  mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = c(0.7, 1, 1))
  lines(density(nonIDs_permTestResults@permuted),
        col = "dodgerblue3",
        lwd = 1.5)
  ablineclip(v = nonIDs_permTestResults@expected,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
  ablineclip(v = nonIDs_permTestResults@observed,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
  ablineclip(v = nonIDs_permTestResults@alpha0.05,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
  text(x = c(textX1,
             nonIDs_permTestResults@expected,
             nonIDs_permTestResults@observed,
             nonIDs_permTestResults@alpha0.05),
       y = c(maxDensityPlus*.95,
             maxDensityPlus,
             maxDensityPlus,
             maxDensityPlus*.95),
       labels = c("Permuted",
                  "Expected",
                  "Observed",
                  expression(alpha*" = 0.05")),
       col = c("dodgerblue",
               "black",
               "forestgreen",
               "red"),
       cex = 0.8)
  dev.off()

  # Determine whether mean orderingFactor values for genes annotated with GO_ID in quantile quantileNo
  # are lower than or higher than those for genes randomly selected from featuresDF_nonIDs and featuresDF_annoGOIDs
  IDs_lessThan_annoGOIDs_Bool <- IDs_mean < ran_annoGOIDs_perm_means
  IDs_moreThan_annoGOIDs_Bool <- IDs_mean > ran_annoGOIDs_perm_means
  
  # Calculate P-values and significance levels
  if(IDs_mean < mean(ran_annoGOIDs_perm_means)) {
    annoGOIDs_pval <- 1 - ( sum(IDs_lessThan_annoGOIDs_Bool) / length(ran_annoGOIDs_perm_means) )
    if(annoGOIDs_pval == 0) {
      annoGOIDs_pval <- minPval
    }
    annoGOIDs_MoreOrLessThanRandom <- "LessThanRandom"
    annoGOIDs_alpha0.05 <- quantile(ran_annoGOIDs_perm_means, probs = 0.05)[[1]]
  } else {
    annoGOIDs_pval <- 1 - ( sum(IDs_moreThan_annoGOIDs_Bool) / length(ran_annoGOIDs_perm_means) )
    if(annoGOIDs_pval == 0) {
      annoGOIDs_pval <- minPval
    }
    annoGOIDs_MoreOrLessThanRandom <- "MoreThanRandom"
    annoGOIDs_alpha0.05 <- quantile(ran_annoGOIDs_perm_means, probs = 0.95)[[1]]
  }
  
  # Create permutation test results object
  annoGOIDs_permTestResults <- new("permTest_popgen_stats",
                                alternative = annoGOIDs_MoreOrLessThanRandom,
                                alpha0.05 = annoGOIDs_alpha0.05,
                                pval = annoGOIDs_pval,
                                observed = IDs_mean,
                                permuted = ran_annoGOIDs_perm_means,
                                expected = mean(ran_annoGOIDs_perm_means),
                                log2obsexp = log2(IDs_mean / mean(ran_annoGOIDs_perm_means)),
                                log2alpha = log2(annoGOIDs_alpha0.05 / mean(ran_annoGOIDs_perm_means)))
  save(annoGOIDs_permTestResults,
       file = paste0(outDir[x],
                     orderingFactor, "_", pop_name[x], "_random_annoGOIDs_permTest_for_", featureNamePlot,
                     "_in_quantile", quantileNo, "_of_", quantiles,
                     "_by_log2_", libName, "_control_in_", region, "_of_",
                     substring(featureName[1][1], first = 1, last = 5), "_in_",
                     paste0(substring(featureName, first = 10, last = 16),
                            collapse = "_"), "_",
                     substring(featureName[1][1], first = 18),
                     "_ann_with_GO_BP_enrichment_GO:", GO_ID, ".RData"))

  # Generate histogram
  pdf(paste0(plotDir[x],
             orderingFactor, "_", pop_name[x], "_random_annoGOIDs_permTest_for_", featureNamePlot,
             "_in_quantile", quantileNo, "_of_", quantiles,
             "_by_log2_", libName, "_control_in_", region, "_of_",
             substring(featureName[1][1], first = 1, last = 5), "_in_",
             paste0(substring(featureName, first = 10, last = 16),
                    collapse = "_"), "_",
             substring(featureName[1][1], first = 18),
             "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_hist.pdf"), 
             height = 4.5, width = 5)
  par(mar = c(3.1, 3.1, 4.1, 1.1),
      mgp = c(1.85, 0.75, 0))
  ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
  #options(scipen = 100)
  # Calculate max density
  maxDensityPlus <- max(density(annoGOIDs_permTestResults@permuted)$y)*1.2
  if(orderingFactor %in% c("TajimaD", "FuLiF", "FuLiD")) {
    if(annoGOIDs_permTestResults@alternative == "MoreThanRandom") {
      xlim <- c(pmin(-1, min(annoGOIDs_permTestResults@permuted)*1.2),
                pmax(annoGOIDs_permTestResults@observed/1.2, annoGOIDs_permTestResults@alpha0.05/1.2))
      textX1 <- quantile(xlim, 0.25)[[1]]
#      textX1 <- min(annoGOIDs_permTestResults@permuted)*1.15
    } else {
      xlim <- c(pmin(-1, annoGOIDs_permTestResults@observed*1.2),
                max(annoGOIDs_permTestResults@permuted)/1.2)
      textX1 <- quantile(xlim, 0.75)[[1]]
#      textX1 <- min(annoGOIDs_permTestResults@permuted)*1.15
    }
  } else {
    if(annoGOIDs_permTestResults@alternative == "MoreThanRandom") {
      xlim <- c(pmin(-1, min(annoGOIDs_permTestResults@permuted)/1.2),
                pmax(annoGOIDs_permTestResults@observed*1.2, annoGOIDs_permTestResults@alpha0.05*1.2))
      textX1 <- quantile(xlim, 0.25)[[1]]
#      textX1 <- min(annoGOIDs_permTestResults@permuted)/1.15
    } else {
      xlim <- c(pmin(-1, annoGOIDs_permTestResults@observed/1.2),
                max(annoGOIDs_permTestResults@permuted)*1.2)
      textX1 <- quantile(xlim, 0.75)[[1]]
#      textX1 <- min(annoGOIDs_permTestResults@permuted)/1.15
    }
  }
  hist(annoGOIDs_permTestResults@permuted,
       freq = FALSE,
       col = "dodgerblue",
       border = NA,
       lwd = 2,
       xlim = c(pretty(xlim)[1],
                pretty(xlim)[length(pretty(xlim))]),
       ylim = c(0,
                maxDensityPlus),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "",
       axes = FALSE)
  axis(side = 2,
       at = pretty(density(annoGOIDs_permTestResults@permuted)$y),
       lwd = 2)
  mtext(side = 2,
        text = "Density",
        line = 1.85)
  axis(side = 1,
       at = pretty(xlim),
       lwd = 2)
  mtext(side = 1,
        text = bquote("Mean" ~ .(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
        line = 1.85)
  titleText <- list(bquote(.(gsub("_", " ", featureNamePlot)) ~
                           "in" ~ .(sub("_\\w+", "", libName)) ~ "Quantile" ~ .(as.character(quantileNo)) ~
                           "(" * .(region) * ") in" ~
                           .(paste0(substring(featureName, first = 10, last = 16),
                                    collapse = " ")) ~ 
                           .(substring(featureName[1][1], first = 18))),
                    bquote(italic("P")*" = "*
                           .(as.character(round(annoGOIDs_permTestResults@pval,
                                                digits = 6)))),
                    bquote("Permutations = "*.(prettyNum(length(annoGOIDs_permTestResults@permuted),
                                                         big.mark = ",",
                                                         trim = T))))
  mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = c(0.7, 1, 1))
  lines(density(annoGOIDs_permTestResults@permuted),
        col = "dodgerblue3",
        lwd = 1.5)
  ablineclip(v = annoGOIDs_permTestResults@expected,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
  ablineclip(v = annoGOIDs_permTestResults@observed,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
  ablineclip(v = annoGOIDs_permTestResults@alpha0.05,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
  text(x = c(textX1,
             annoGOIDs_permTestResults@expected,
             annoGOIDs_permTestResults@observed,
             annoGOIDs_permTestResults@alpha0.05),
       y = c(maxDensityPlus*.95,
             maxDensityPlus,
             maxDensityPlus,
             maxDensityPlus*.95),
       labels = c("Permuted",
                  "Expected",
                  "Observed",
                  expression(alpha*" = 0.05")),
       col = c("dodgerblue",
               "black",
               "forestgreen",
               "red"),
       cex = 0.8)
  dev.off()
}
#}, mc.cores = length(pop_name), mc.preschedule = F)
