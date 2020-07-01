#!/applications/R/R-3.5.0/bin/Rscript

# Compare population genetics statistics for GO-term-annotated genes in a given ASY1 quantile
# with 1) those for GO-term-annotated genes in another given ASY1 quantile (using LSD, t-tests, and Yuen t-tests/MWW tests), and
#      2) those for randomSets sets of randomly selected genes in another given ASY1 quantile and not encoding NLRs (using permutation tests)

# Usage:
# /applications/R/R-3.5.0/bin/Rscript quantile_genes_response_to_cold_GO_genes_popgen_stats_permTest.R ASY1_CS_Rep1_ChIP ASY1_CS 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' 'Cold_response_genes' '0009409' genes 1 3 4 TajimaD_all "Tajima's D" 10000 0.0001 '%4.0f' '4.1,3.5'

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#featureNamePlot <- "Cold_response_genes"
#GO_ID <- "0009409"
#region <- "genes"
#quantileNo <- 1
#firstLowerQuantile <- 3
#quantiles <- 4
#orderingFactor <- "TajimaD_all"
#orderingFactor <- "RozasR2_all"
#orderingFactor <- "CLR_all"
#orderingFactorName <- bquote("Tajima's" ~ italic("D"))
#orderingFactorName <- bquote("Rozas'" ~ italic("R")[2])
#orderingFactorName <- unlist(strsplit("Tajima's D", split = " "))
#orderingFactorName <- unlist(strsplit("Rozas' R 2", split = " "))
#randomSets <- 10000
#minPval <- 0.0001
## For Tajima's D
#yDec <- "%4.0f"
#xAnn <- as.numeric(unlist(strsplit("4.1,3.5", split = ",")))
## For Rozas' R2
#yDec <- "%1.2f"
#xAnn <- as.numeric(unlist(strsplit("0.26,0.23", split = ",")))
## For CLR
#yDec <- "%4.0f"
#xAnn <- as.numeric(unlist(strsplit("350,325", split = ",")))

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
featureNamePlot <- args[4]
GO_ID <- args[5]
region <- args[6]
quantileNo <- as.numeric(args[7])
firstLowerQuantile <- as.numeric(args[8])
quantiles <- as.numeric(args[9])
orderingFactor <- args[10]
orderingFactorName <- unlist(strsplit(args[11], split = " "))
if(grepl("Tajima", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ italic(.(orderingFactorName[2])))
} else if(grepl("Rozas' Z", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ italic(.(orderingFactorName[2])))
} else if(grepl("Rozas' R", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ italic(.(orderingFactorName[2]))[.(as.numeric(orderingFactorName[3]))])
} else if(grepl("Diversity", paste(orderingFactorName, collapse = " "))) {
  orderingFactorName <- bquote(.(orderingFactorName[1]) ~ "(" * .(as.symbol(orderingFactorName[2])) * ")")
} else {
  orderingFactorName <- paste(orderingFactorName, collapse = " ")
}
randomSets <- as.numeric(args[12])
minPval <- as.numeric(args[13])
yDec <- as.character(args[14])
xAnn <- as.numeric(unlist(strsplit(args[15], split = ",")))

library(parallel)
library(plotrix)
#library(tidyr)
library(dplyr)
library(WRS2)
library(PairedData)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

pop_name <- c(
              "Africa",
              "MiddleEast",
              "Asia",
              "FormerSU",
              "EasternEurope",
              "WesternEurope",
              "NorthAmerica",
              "CentralAmerica",
              "SouthAmerica",
              "Oceania"
             )

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

if(libName %in% c("cMMb", "HudsonRM_all")) {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName), "/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", region, "/")
}
outDir <- sapply(seq_along(pop_name), function(x) {
  paste0(outDir, pop_name[x], "/")
})
plotDir <- paste0(outDir, "plots/")

# Define quantile colours
quantileColours <- c("red", "navy")
makeTransparent <- function(thisColour, alpha = 180)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
quantileColours <- makeTransparent(quantileColours)

# Disable scientific notation (e.g., 0.0001 rather than 1e-04)
options(scipen = 100)

# Load functional annotation in order to extract "response to cold" genes
anno <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/RamirezGonzalez_2018_Science_GO_anno/",
                          "RamirezGonzalez_2018_iwgsc_refseqv1.0_OntologiesForGenes_FunctionalAnnotation_HCgenes_in_Agenome_Bgenome_Dgenome_genomewide_GO_IDs_no_chrUn.tsv"),
                   sep = "\t", stringsAsFactors = F)
colnames(anno) <- c("featureID", "GO")
response_to_cold_indices <- which(grepl(pattern = GO_ID, x = anno$GO))
response_to_cold_indices <- sort(unique(c(response_to_cold_indices)))
# Retain only "response to cold" genes
# Get "response to cold" gene IDs and their row indices in features
IDs <- anno[response_to_cold_indices,]$featureID

# Load table of features grouped into quantiles
# by decreasing log2(libName/control) in region
#mclapply(seq_along(pop_name), function(x) {
estimates_allpops <- data.frame()
IDsDF_annoGOIDsDF_stat_allpops <- data.frame() 
for(x in seq_along(pop_name)) {
  print(pop_name[x])
  if(libName %in% c("cMMb", "HudsonRM_all")) {
  featuresDF <- read.table(paste0(outDir[x], "features_", quantiles, "quantiles",
                                  "_by_", sub("_\\w+", "", libName), "_of_",
                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
                                  paste0(substring(featureName, first = 10, last = 16),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                           header = T, sep = "\t", stringsAsFactors = F)
  } else {
  featuresDF <- read.table(paste0(outDir[x], "features_", quantiles, "quantiles",
                                  "_by_", sub("_\\w+", "", libName), "_in_",
                                  region, "_of_",
                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
                                  paste0(substring(featureName, first = 10, last = 16),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                           header = T, sep = "\t", stringsAsFactors = F)
  }

  featuresDF$featureID <- sub("\\.\\d+", "", featuresDF$featureID)
  # Retain rows that contain non-NA values of orderingFactor
  featuresDF <- featuresDF[!is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor)]),]
  
  # Subset genes (GO-term-annotated genes in Quantile 1 [IDs],
  # GO-term-annotated genes in, e.g., Quantile 3 and Quantile 4 (annoGOIDs), and
  # genes not annotated with given GO term and in, e.g., Quantile 3 and Quantile 4 [nonIDs])
  # Below is commented out in favour of comparison of a more exhaustive set of genes
#  featuresDF_IDsDF <- featuresDF[featuresDF$featureID %in% IDs,]
  featuresDF_IDsDF <- featuresDF[featuresDF$featureID %in% IDs &
                                 featuresDF$quantile == paste0("Quantile ", quantileNo),]
#  featuresDF_nonIDsDF <- featuresDF[!(featuresDF$featureID %in% IDs),]
#  featuresDF_nonIDsDF <- featuresDF[featuresDF$quantile != paste0("Quantile ", quantileNo),]
#  featuresDF_nonIDsDF <- featuresDF[featuresDF$quantile == paste0("Quantile ", quantiles),]
  featuresDF_nonIDsDF <- featuresDF[!(featuresDF$featureID %in% IDs) &
                                    featuresDF$quantile %in% c(paste0("Quantile ", firstLowerQuantile), paste0("Quantile ", quantiles)),]
  featuresDF_annoGOIDsDF <- featuresDF[featuresDF$featureID %in% IDs &
                                       featuresDF$quantile %in% c(paste0("Quantile ", firstLowerQuantile), paste0("Quantile ", quantiles)),]
  featuresDF_annoGOIDsDF$quantile <- paste0("Quantiles ", firstLowerQuantile, " & ", quantiles)
 
  ### NOTE THAT PRE-TRIMMING THE DATA MAY MAKE YUEN T-TESTS BELOW INVALID
  ### APPLIED HERE AS A TEST DUE TO RozasR2 OUTLIERS;
  ### REMOVING THEM MAKES PLOTS EASIER TO INTERPRET
  if(grepl("RozasR2", orderingFactor)) {
   featuresDF_IDsDF <-  featuresDF_IDsDF[featuresDF_IDsDF[,which(colnames(featuresDF_IDsDF) ==
                                                                 orderingFactor)] < 1,]
   featuresDF_annoGOIDsDF <-  featuresDF_annoGOIDsDF[featuresDF_annoGOIDsDF[,which(colnames(featuresDF_annoGOIDsDF) ==
                                                                                   orderingFactor)] < 1,]
   featuresDF_nonIDsDF <-  featuresDF_nonIDsDF[featuresDF_nonIDsDF[,which(colnames(featuresDF_nonIDsDF) ==
                                                                          orderingFactor)] < 1,]
  }
  ### NOTE THAT PRE-TRIMMING THE DATA MAY MAKE YUEN T-TESTS BELOW INVALID

  featuresDF_IDs <- featuresDF_IDsDF$featureID
  featuresDF_nonIDs <- featuresDF_nonIDsDF$featureID
  featuresDF_annoGOIDs <- featuresDF_annoGOIDsDF$featureID

  # Combine featuresDF_IDsDF and featuresDF_annoGOIDsDF to enable calculation of LSDs
  IDsDF_annoGOIDsDF <- rbind(featuresDF_IDsDF, featuresDF_annoGOIDsDF)
  IDsDF_annoGOIDsDF_stat <- data.frame(population = pop_name_plot[x],
                                       stat = IDsDF_annoGOIDsDF[,which(colnames(IDsDF_annoGOIDsDF) == orderingFactor)],
                                       quantile = IDsDF_annoGOIDsDF$quantile)
  colnames(IDsDF_annoGOIDsDF_stat) <- c("population", orderingFactor, "quantile")
  IDsDF_annoGOIDsDF_stat_allpops <- rbind(IDsDF_annoGOIDsDF_stat_allpops, IDsDF_annoGOIDsDF_stat)

  # Linear model
#  lm1 <- lm(TajimaD ~ quantile, data = IDsDF_annoGOIDsDF)
  lm1 <- lm(IDsDF_annoGOIDsDF[,which(colnames(IDsDF_annoGOIDsDF) == orderingFactor)] ~
            IDsDF_annoGOIDsDF$quantile)
  # Create a dataframe containing means, SDs, SEMs or interval bounds
  estimates <- expand.grid(population = pop_name_plot[x],
                           quantile = unique(IDsDF_annoGOIDsDF$quantile)) 
  # Add the mean orderingFactor values to the dataframe
  estimates$mean <- c(mean(IDsDF_annoGOIDsDF[IDsDF_annoGOIDsDF$quantile == paste0("Quantile ", quantileNo),
                                             which(colnames(IDsDF_annoGOIDsDF) == orderingFactor)]),
                      mean(IDsDF_annoGOIDsDF[IDsDF_annoGOIDsDF$quantile == paste0("Quantiles ", firstLowerQuantile, " & ", quantiles),
                                             which(colnames(IDsDF_annoGOIDsDF) == orderingFactor)]))
  # Add the standard error of the difference between means to the estimates dataframe
  estimates$sed <- summary(lm1)$coefficients[2,2]
  # Get gene-grouping-specific 95% critical value of t for the gene-grouping-specific df
  alpha <- 0.05
  tQuantile1 <- qt(p = 1-(alpha/2),
                   df = dim(IDsDF_annoGOIDsDF[IDsDF_annoGOIDsDF$quantile == paste0("Quantile ", quantileNo),])[1] - 1)
  tQuantile4 <- qt(p = 1-(alpha/2),
                   df = dim(IDsDF_annoGOIDsDF[IDsDF_annoGOIDsDF$quantile == paste0("Quantiles ", firstLowerQuantile, " & ", quantiles),])[1] - 1) 
  estimates$lsd <- c(tQuantile1*estimates$sed[1], tQuantile4*estimates$sed[2])

  if(estimates$mean[1] < estimates$mean[2]) {
    altHyp <- "less"
  } else if(estimates$mean[1] > estimates$mean[2]) {
    altHyp <- "greater"
  } else {
    altHyp <- "two.sided"
  }

  # Evaluate differences between gene-quantile orderingFactor values
  Utest <- wilcox.test(x = featuresDF_IDsDF[,which(colnames(featuresDF_IDsDF) ==
                                                   orderingFactor)],
                       y = featuresDF_annoGOIDsDF[,which(colnames(featuresDF_annoGOIDsDF) ==                  
                                                         orderingFactor)],
                       alternative = altHyp)
  UtestPval <- Utest$p.value
  UtestPvalChar <- if(UtestPval < 0.0001) {
                     "< 0.0001"
                   } else {
                     paste0("= ", as.character(round(UtestPval, digits = 4)))
                   }

  ttest <- stats::t.test(x = featuresDF_IDsDF[,which(colnames(featuresDF_IDsDF) ==
                                                     orderingFactor)],
                         y = featuresDF_annoGOIDsDF[,which(colnames(featuresDF_annoGOIDsDF) ==                  
                                                           orderingFactor)],
                         alternative = altHyp)
  ttestPval <- ttest$p.value
  ttestPvalChar <- if(ttestPval < 0.0001) {
                     "< 0.0001"
                   } else {
                     paste0("= ", as.character(round(ttestPval, digits = 4)))
                   }

  trim <- 0.1
##  yuenbttest <- yuenbt(TajimaD ~ quantile, data = IDsDF_annoGOIDsDF,
##                       tr = trim, nboot = 10000, side = T)
#  yuenbttest <- yuenbt(IDsDF_annoGOIDsDF[,which(colnames(IDsDF_annoGOIDsDF) == orderingFactor)] ~
#                       IDsDF_annoGOIDsDF$quantile, 
#                       tr = trim, nboot = 10000, side = T)
#  yuenbttestPval <- yuenbttest$p.value
#  yuenbttestPval <- yuenbttestPval/2
#  yuenbttestPvalChar <- if(yuenbttestPval < 0.0001) {
#                          "< 0.0001"
#                        } else {
#                          paste0("= ", as.character(round(yuenbttestPval, digits = 4)))
#                        }
  # Bootstrapped version above does not provide for one-tailed tests, so using yuen.t.test instead
  # Alternatively, could divide yuenbttestPval by 2
  yuenttest <- yuen.t.test(x = featuresDF_IDsDF[,which(colnames(featuresDF_IDsDF) ==
                                                       orderingFactor)],
                           y = featuresDF_annoGOIDsDF[,which(colnames(featuresDF_annoGOIDsDF) ==
                                                       orderingFactor)],
                           tr = trim, alternative = altHyp)
  yuenttestPval <- yuenttest$p.value
  yuenttestPvalChar <- if(yuenttestPval < 0.0001) {
                         "< 0.0001"
                       } else {
                         paste0("= ", as.character(round(yuenttestPval, digits = 4)))
                       }

  estimates$UtestPval <- UtestPvalChar
  estimates$ttestPval <- ttestPvalChar
  estimates$yuenttestPval <- yuenttestPvalChar
#  estimates$yuenbttestPval <- yuenbttestPvalChar

  # Combine estimates for each population into one dataframe 
  estimates_allpops <- rbind(estimates_allpops, estimates)

  # Plot orderingFactor means and LSDs for IDs vs annoGOIDs
  popgen_stats_meanLSDs <- function(dataFrame,
                                    parameterLab,
                                    featureGroup,
                                    featureNamePlot) {
    ggplot(data = dataFrame,
           mapping = aes(x = get(featureGroup),
                         y = mean,
                         colour = get(featureGroup))) +
    labs(colour = "") +
    geom_point(shape = "-", size = 18, position = position_dodge(width = 0.2), colour = "black") +
    geom_errorbar(mapping = aes(ymin = mean-(lsd/2),
                                ymax = mean+(lsd/2)),
                  width = 0.5, size = 1, position = position_dodge(width = 0.2), colour = "black") +
    geom_beeswarm(data = IDsDF_annoGOIDsDF,
                  mapping = aes(x = get(featureGroup),
                                y = get(orderingFactor),
                                colour = get(featureGroup)),
                  priority = "ascending",
                  cex = 1,
                  size = 1) +
    scale_colour_manual(values = quantileColours) +
    scale_y_continuous(
#                       limits = c(summary_stats_min, summary_stats_max),
                       labels = function(x) sprintf(yDec, x)) +
#    scale_x_discrete(breaks = as.vector(dataFrame$quantile),
#                     labels = as.vector(dataFrame$quantile)) +
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
          plot.margin = unit(c(0.8,1.2,0.1,0.3),"cm"),
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote(atop(.(featureNamePlot),
                        atop(italic("t") * "-test" ~ italic("P") ~ .(ttestPvalChar),
                             "Yuen's test (10% trimmed mean)"  ~ italic("P") ~ .(yuenttestPvalChar)))))
#                             "MWW test"  ~ italic("P") ~ .(UtestPvalChar)))))
  }
  
  ggObjGA_feature_mean <- popgen_stats_meanLSDs(dataFrame = estimates,
                                                parameterLab = bquote(.(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
                                                featureGroup = "quantile",
                                                featureNamePlot = gsub("_", " ", featureNamePlot)
                                               )
  
  if(libName %in% c("cMMb", "HudsonRM_all")) {
  ggsave(paste0(plotDir[x],
                orderingFactor, "_", pop_name[x], "_IDs_v_annoGOIDs_for_", gsub(" ", "_", featureNamePlot),
                "_in_quantile", quantileNo, "_of_", quantiles,
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18),
                "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_meanLSD_v010720.pdf"),
         plot = ggObjGA_feature_mean,
         height = 8, width = 7)
  } else {
  ggsave(paste0(plotDir[x],
                orderingFactor, "_", pop_name[x], "_IDs_v_annoGOIDs_for_", gsub(" ", "_", featureNamePlot),
                "_in_quantile", quantileNo, "_of_", quantiles,
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18),
                "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_meanLSD_v010720.pdf"),
         plot = ggObjGA_feature_mean,
         height = 8, width = 7)
  }

#  # Calculate mean orderingFactor value for genes annotated with GO_ID in quantile quantileNo
#  IDs_mean <- mean(featuresDF_IDsDF[,which(colnames(featuresDF_IDsDF) ==
#                                           orderingFactor)],
#                   na.rm = T)
#  
#  # Function to randomly select n feature IDs not present in IDs (from featuresDF_nonIDs or featuresDF_annoGOIDs)
#  ran_nonIDs_select <- function(nonIDsChr, n) {
#    sample(x = nonIDsChr,
#           size = n,
#           replace = FALSE)
#  }
#  
#  # Define seed so that random selections are reproducible
#  set.seed(453838430)
#  
#  # Apply ran_nonIDs_select() function on a per-chromosome basis to generate randomSets random sets
#  # of IDs from featuresDF_nonIDs for which to calculate randomSets mean orderingFactor values
#  ran_nonIDs_perm_means <- unlist(mclapply(1:randomSets, function(y) {
#    ran_nonIDs <- NULL
#    for(i in 1:length(unique(featuresDF_IDsDF$seqnames))) {
#      IDsChr <- featuresDF_IDs[grepl(paste0("TraesCS",
#                                            sub("chr", "", unique(featuresDF_IDsDF$seqnames))[i]),
#                                     featuresDF_IDs)]
#      nonIDsChr <- featuresDF_nonIDs[grepl(paste0("TraesCS",
#                                                  sub("chr", "", unique(featuresDF_IDsDF$seqnames))[i]),
#                                           featuresDF_nonIDs)]
#      ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
#                                         n = length(IDsChr))
#      ran_nonIDs <- c(ran_nonIDs, ran_nonIDsChr)
#    }
#    return(mean(featuresDF_nonIDsDF[featuresDF_nonIDsDF$featureID %in%
#                                    ran_nonIDs,][,which(colnames(featuresDF_nonIDsDF) ==
#                                                        orderingFactor)],
#                na.rm = T))
#  }, mc.cores = detectCores(), mc.preschedule = T))
#  
#  # Set class for permutation test results object
#  setClass("permTest_popgen_stats",
#           representation(alternative = "character",
#                          alpha0.05 = "numeric",
#                          pval = "numeric",
#                          observed = "numeric",
#                          permuted = "numeric",
#                          expected = "numeric",
#                          log2obsexp = "numeric",
#                          log2alpha = "numeric",
#                          GOgenesInQuantile = "numeric",
#                          GOgenesNotInQuantile = "numeric",
#                          genesNotInQuantile = "numeric"))
#  
#  # Disable scientific notation (e.g., 0.0001 rather than 1e-04)
#  options(scipen = 100)
#  
#  # Determine whether mean orderingFactor values for genes annotated with GO_ID in quantile quantileNo
#  # are lower than or higher than those for genes randomly selected from featuresDF_nonIDs and featuresDF_annoGOIDs
#  IDs_lessThan_nonIDs_Bool <- IDs_mean < ran_nonIDs_perm_means
#  IDs_moreThan_nonIDs_Bool <- IDs_mean > ran_nonIDs_perm_means
#  
#  # Calculate P-values and significance levels
#  if(IDs_mean < mean(ran_nonIDs_perm_means)) {
#    nonIDs_pval <- 1 - ( sum(IDs_lessThan_nonIDs_Bool) / length(ran_nonIDs_perm_means) )
#    if(nonIDs_pval == 0) {
#      nonIDs_pval <- minPval
#    }
#    nonIDs_MoreOrLessThanRandom <- "LessThanRandom"
#    nonIDs_alpha0.05 <- quantile(ran_nonIDs_perm_means, probs = 0.05)[[1]]
#  } else {
#    nonIDs_pval <- 1 - ( sum(IDs_moreThan_nonIDs_Bool) / length(ran_nonIDs_perm_means) )
#    if(nonIDs_pval == 0) {
#      nonIDs_pval <- minPval
#    }
#    nonIDs_MoreOrLessThanRandom <- "MoreThanRandom"
#    nonIDs_alpha0.05 <- quantile(ran_nonIDs_perm_means, probs = 0.95)[[1]]
#  }
#  
#  # Create permutation test results object
#  nonIDs_permTestResults <- new("permTest_popgen_stats",
#                                alternative = nonIDs_MoreOrLessThanRandom,
#                                alpha0.05 = nonIDs_alpha0.05,
#                                pval = nonIDs_pval,
#                                observed = IDs_mean,
#                                permuted = ran_nonIDs_perm_means,
#                                expected = mean(ran_nonIDs_perm_means),
#                                log2obsexp = log2(IDs_mean / mean(ran_nonIDs_perm_means)),
#                                log2alpha = log2(nonIDs_alpha0.05 / mean(ran_nonIDs_perm_means)),
#                                GOgenesInQuantile = length(featuresDF_IDs),
#                                GOgenesNotInQuantile = length(featuresDF_annoGOIDs),
#                                genesNotInQuantile = length(featuresDF_nonIDs))
#  if(libName %in% c("cMMb", "HudsonRM_all")) {
#  save(nonIDs_permTestResults,
#       file = paste0(outDir[x],
#                     orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#                     "_in_quantile", quantileNo, "_of_", quantiles,
#                     "_by_", libName, "_of_",
#                     substring(featureName[1][1], first = 1, last = 5), "_in_",
#                     paste0(substring(featureName, first = 10, last = 16),
#                            collapse = "_"), "_",
#                     substring(featureName[1][1], first = 18),
#                     "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_v010720.RData"))
#  } else {
#  save(nonIDs_permTestResults,
#       file = paste0(outDir[x],
#                     orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#                     "_in_quantile", quantileNo, "_of_", quantiles,
#                     "_by_log2_", libName, "_control_in_", region, "_of_",
#                     substring(featureName[1][1], first = 1, last = 5), "_in_",
#                     paste0(substring(featureName, first = 10, last = 16),
#                            collapse = "_"), "_",
#                     substring(featureName[1][1], first = 18),
#                     "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_v010720.RData"))
#  }
#
#  # Generate histogram
#  if(libName %in% c("cMMb", "HudsonRM_all")) {
#  pdf(paste0(plotDir[x],
#             orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#             "_in_quantile", quantileNo, "_of_", quantiles,
#             "_by_", libName, "_of_",
#             substring(featureName[1][1], first = 1, last = 5), "_in_",
#             paste0(substring(featureName, first = 10, last = 16),
#                    collapse = "_"), "_",
#             substring(featureName[1][1], first = 18),
#             "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_hist_v010720.pdf"),
#             height = 4.5, width = 5)
#  } else {
#  pdf(paste0(plotDir[x],
#             orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#             "_in_quantile", quantileNo, "_of_", quantiles,
#             "_by_log2_", libName, "_control_in_", region, "_of_",
#             substring(featureName[1][1], first = 1, last = 5), "_in_",
#             paste0(substring(featureName, first = 10, last = 16),
#                    collapse = "_"), "_",
#             substring(featureName[1][1], first = 18),
#             "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_hist_v010720.pdf"),
#             height = 4.5, width = 5)
#  }
#  par(mar = c(3.1, 3.1, 4.1, 1.1),
#      mgp = c(1.85, 0.75, 0))
#  ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
#  #options(scipen = 100)
#  # Calculate max density
#  maxDensityPlus <- max(density(nonIDs_permTestResults@permuted)$y)*1.2
#  if(orderingFactor %in% c("TajimaD_all", "TajimaD_syn", "TajimaD_nonsysn",
#                           "FuLiF_all", "FuLiF_syn", "FuLiF_nonsyn",
#                           "FuLiD_all", "FuLiD_syn", "FuLiD_nonsyn")) {
#    if(nonIDs_permTestResults@alternative == "MoreThanRandom") {
#      xlim <- c(pmin(-1, min(nonIDs_permTestResults@permuted)*1.2),
#                pmax(nonIDs_permTestResults@observed/1.2, nonIDs_permTestResults@alpha0.05/1.2))
#      textX1 <- quantile(xlim, 0.25)[[1]]
##      textX1 <- min(nonIDs_permTestResults@permuted)*1.15
#    } else {
#      xlim <- c(pmin(-1, nonIDs_permTestResults@observed*1.2),
#                max(nonIDs_permTestResults@permuted)/1.2)
#      textX1 <- quantile(xlim, 0.75)[[1]]
##      textX1 <- min(nonIDs_permTestResults@permuted)*1.15
#    }
#  } else {
#    if(nonIDs_permTestResults@alternative == "MoreThanRandom") {
#      xlim <- c(pmin(0, min(nonIDs_permTestResults@permuted)/1.2),
#                pmax(nonIDs_permTestResults@observed*1.2, nonIDs_permTestResults@alpha0.05*1.2))
#      textX1 <- quantile(xlim, 0.25)[[1]]
##      textX1 <- min(nonIDs_permTestResults@permuted)/1.15
#    } else {
#      xlim <- c(pmin(0, nonIDs_permTestResults@observed/1.2),
#                max(nonIDs_permTestResults@permuted)*1.2)
#      textX1 <- quantile(xlim, 0.75)[[1]]
##      textX1 <- min(nonIDs_permTestResults@permuted)/1.15
#    }
#  }
#  hist(nonIDs_permTestResults@permuted,
#       breaks = 50,
#       freq = FALSE,
#       col = "dodgerblue",
#       border = NA,
#       lwd = 2,
#       xlim = c(pretty(xlim)[1],
#                pretty(xlim)[length(pretty(xlim))]),
#       ylim = c(0,
#                maxDensityPlus),
#       xaxt = "n", yaxt = "n",
#       xlab = "", ylab = "", main = "",
#       axes = FALSE)
#  axis(side = 2,
#       at = pretty(density(nonIDs_permTestResults@permuted)$y),
#       lwd = 2)
#  mtext(side = 2,
#        text = "Density",
#        line = 1.85)
#  axis(side = 1,
#       at = pretty(xlim),
#       lwd = 2)
#  mtext(side = 1,
#        text = bquote("Mean" ~ .(orderingFactorName) ~ "(" * .(pop_name_plot[x]) * ")"),
#        line = 1.85)
#  titleText <- list(bquote(.(gsub("_", " ", featureNamePlot)) ~
#                           "in" ~ .(sub("_\\w+", "", libName)) ~ "Quantile" ~ .(as.character(quantileNo)) ~
#                           "(" * .(region) * ") in" ~
#                           .(paste0(substring(featureName, first = 10, last = 16),
#                                    collapse = " ")) ~ 
#                           .(substring(featureName[1][1], first = 18))),
#                    bquote(italic("P")*" = "*
#                           .(as.character(round(nonIDs_permTestResults@pval,
#                                                digits = 6)))),
#                    bquote("Permutations = "*.(prettyNum(length(nonIDs_permTestResults@permuted),
#                                                         big.mark = ",",
#                                                         trim = T))))
#  mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = c(0.7, 1, 1))
#  lines(density(nonIDs_permTestResults@permuted),
#        col = "dodgerblue3",
#        lwd = 1.5)
#  ablineclip(v = nonIDs_permTestResults@expected,
#             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
#  ablineclip(v = nonIDs_permTestResults@observed,
#             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
#  ablineclip(v = nonIDs_permTestResults@alpha0.05,
#             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
#  text(x = c(textX1,
#             nonIDs_permTestResults@expected,
#             nonIDs_permTestResults@observed,
#             nonIDs_permTestResults@alpha0.05),
#       y = c(maxDensityPlus*.95,
#             maxDensityPlus,
#             maxDensityPlus,
#             maxDensityPlus*.95),
#       labels = c("Permuted",
#                  "Expected",
#                  "Observed",
#                  expression(alpha*" = 0.05")),
#       col = c("dodgerblue",
#               "black",
#               "forestgreen",
#               "red"),
#       cex = 0.8)
#  dev.off()
}

# Plot orderingFactor means and LSDs for IDs vs annoGOIDs for all populations
popgen_stats_meanLSDs <- function(dataFrame1,
                                  dataFrame2,
                                  parameterLab,
                                  populationGroup,
                                  featureGroup,
                                  featureNamePlot) {
  ggplot(data = dataFrame1,
         mapping = aes(x = get(populationGroup),
                       y = mean,
                       group = get(featureGroup))) +
  labs(colour = "") +
  geom_point(shape = "-", size = 18, position =  position_dodge(width = 0.6)) +
  geom_errorbar(mapping = aes(ymin = mean-(lsd/2),
                              ymax = mean+(lsd/2)),
                width = 0.5, size = 2, position =  position_dodge(width = 0.6)) +
  geom_beeswarm(data = dataFrame2,
                mapping = aes(x = get(populationGroup),
                              y = get(orderingFactor),
                              colour = get(featureGroup)),
                priority = "ascending",
                groupOnX = T,
                dodge.width = 0.6,
                cex = 0.2,
                size = 1) +
  scale_colour_manual(values = quantileColours) +
  scale_y_continuous(labels = function(x) sprintf(yDec, x)) +
#  scale_x_discrete(position = "bottom",
#                   breaks = levels(dataFrame1$population),
#                   labels = levels(dataFrame1$population)) +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "top",
                             label.theme = element_text(size = 22, hjust = 0, vjust = 0.5, angle = 90),
                             nrow = 1,
                             byrow = TRUE),
         colour = guide_legend(override.aes = list(size = 10))) +
  labs(x = "",
       y = parameterLab) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.y = element_line(size = 2.0, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(size = 22, colour = "black", family = "Luxi Mono"),
        axis.text.x = element_text(size = 22, colour = "black", hjust = 0.5, vjust = 1.0, angle = 0),
        axis.title = element_text(size = 26, colour = "black"),
        legend.text = element_text(size = 26),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.8,1.2,0.1,0.3),"cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot))) +
  annotate(geom = "text",
           size = 8,
           x = seq_along(pop_name),
           y = xAnn[1], angle = 0,
           parse = T,
           label = lapply(sapply(seq_along(pop_name)-1, function(w) { (w+1)+w }),
                     function(x) { bquote(italic("t") * "-test" ~ italic("P") ~ .(dataFrame1$ttestPval[x])) }) ) +
  annotate(geom = "text",
           size = 8,
           x = seq_along(pop_name),
           y = xAnn[2], angle = 0,
           parse = T,
           label = lapply(sapply(seq_along(pop_name)-1, function(w) { (w+1)+w }),
                     function(x) { bquote("Yuen" ~ italic("P") ~ .(dataFrame1$yuenttestPval[x])) }) )
#           label = lapply(sapply(seq_along(pop_name)-1, function(w) { (w+1)+w }),
#                     function(x) { bquote("MWW" ~ italic("P") ~ .(dataFrame1$UtestPval[x])) }) )
}
ggObjGA_feature_mean <- popgen_stats_meanLSDs(dataFrame1 = estimates_allpops,
                                              dataFrame2 = IDsDF_annoGOIDsDF_stat_allpops,
                                              parameterLab = bquote(.(orderingFactorName)),
                                              populationGroup = "population",
                                              featureGroup = "quantile",
                                              featureNamePlot = gsub("_", " ", featureNamePlot))
if(libName %in% c("cMMb", "HudsonRM_all")) {
ggsave(paste0(sub("\\w+\\/$", "", outDir[1]),
              orderingFactor, "_allpops_IDs_v_annoGOIDs_for_", gsub(" ", "_", featureNamePlot),
              "_in_quantile", quantileNo, "_of_", quantiles,
              "_by_", libName, "_of_",
              substring(featureName[1][1], first = 1, last = 5), "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18),
              "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_meanLSD_v010720.pdf"),
       plot = ggObjGA_feature_mean,
       height = 9, width = 35)
} else {
ggsave(paste0(sub("\\w+\\/$", "", outDir[1]),
              orderingFactor, "_allpops_IDs_v_annoGOIDs_for_", gsub(" ", "_", featureNamePlot),
              "_in_quantile", quantileNo, "_of_", quantiles,
              "_by_log2_", libName, "_control_in_", region, "_of_",
              substring(featureName[1][1], first = 1, last = 5), "_in_",
              paste0(substring(featureName, first = 10, last = 16),
                     collapse = "_"), "_",
              substring(featureName[1][1], first = 18),
              "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_meanLSD_v010720.pdf"),
       plot = ggObjGA_feature_mean,
       height = 9, width = 35)
}

## Plot bar graph summarising permutation test results
#pt_list <- list()
#for(x in seq_along(pop_name)) {
#  if(libName %in% c("cMMb", "HudsonRM_all")) {
#  load(paste0(outDir[x],
#              orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#              "_in_quantile", quantileNo, "_of_", quantiles,
#              "_by_", libName, "_of_",
#              substring(featureName[1][1], first = 1, last = 5), "_in_",
#              paste0(substring(featureName, first = 10, last = 16),
#                     collapse = "_"), "_",
#              substring(featureName[1][1], first = 18),
#              "_v010720.RData"))
#  pt_list <- c(pt_list, nonIDs_permTestResults)
#  } else {
#  load(paste0(outDir[x],
#              orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#              "_in_quantile", quantileNo, "_of_", quantiles,
#              "_by_log2_", libName, "_control_in_", region, "_of_",
#              substring(featureName[1][1], first = 1, last = 5), "_in_",
#              paste0(substring(featureName, first = 10, last = 16),
#                     collapse = "_"), "_",
#              substring(featureName[1][1], first = 18),
#              "_v010720.RData"))
#  pt_list <- c(pt_list, nonIDs_permTestResults)
# }
#}
#bargraph_df <- data.frame(population = pop_name_plot,
#                          log2ObsExp = sapply(seq_along(pt_list), function(x) { pt_list[[x]]@log2obsexp }),
#                          log2alpha0.05 = sapply(seq_along(pt_list), function(x) { pt_list[[x]]@log2alpha }))
#bargraph_df$population <- factor(bargraph_df$population,
#                                 levels = pop_name_plot)
#bp <- ggplot(data = bargraph_df,
#             mapping = aes(x = population,
#                           y = log2ObsExp,
#                           fill = " ")) +
#  geom_bar(stat = "identity",
#           position = position_dodge()) +
#  scale_fill_manual(name = "",
#                    values = c("dodgerblue3"),
#                    labels = " ") +
#  geom_point(mapping = aes(x = population,
#                           y = log2alpha0.05),
#             position = position_dodge(0.9),
#             shape = "-", colour  = "grey80", size = 20) +
#  labs(y = bquote("Log"[2]*"(observed/expected)" ~ .(orderingFactorName))) +
##  scale_y_continuous(limits = c(-1.5, 1.5)) +
#  scale_x_discrete(position = "top") +
#  guides(fill = guide_legend(direction = "horizontal",
#                             label.position = "top",
#                             label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
#                             nrow = 1,
#                             byrow = TRUE)) +
#  theme_bw() +
#  theme(axis.line.y = element_line(size = 1, colour = "black"),
#        axis.ticks.y = element_line(size = 1, colour = "black"),
#        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
#        axis.title.y = element_text(size = 20, colour = "black"),
#        axis.ticks.x = element_blank(),
#        axis.text.x = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, angle = 90),
#        axis.title.x = element_blank(),
#        panel.grid = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        legend.position = "none",
#        #legend.position = c(0.05, 0.30),
#        legend.background = element_rect(fill = "transparent"),
#        legend.key = element_rect(colour = "transparent",
#                                  fill = "transparent"),
#        plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
#        plot.title = element_text(size = 16, colour = "black", hjust = 0.5)) +
#  ggtitle(bquote(.(gsub("_", " ", featureNamePlot)) ~
#                 "in" ~ .(sub("_\\w+", "", libName)) ~ "Quantile" ~ .(as.character(quantileNo)) ~
#                 "(" * .(region) * ") in" ~
#                 .(paste0(substring(featureName, first = 10, last = 16),
#                          collapse = " ")) ~
#                 .(substring(featureName[1][1], first = 18)) ~
#                 "(" * .(prettyNum(as.character(randomSets),
#                                   big.mark = ",", trim = "T")) ~ "permutations)"))
#if(libName %in% c("cMMb", "HudsonRM_all")) {
#ggsave(paste0(sub("\\w+\\/$", "", outDir[1]),
#              orderingFactor, "_bargraph_allpops_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#              "_in_quantile", quantileNo, "_of_", quantiles,
#              "_by_", libName, "_of_",
#              substring(featureName[1][1], first = 1, last = 5), "_in_",
#              paste0(substring(featureName, first = 10, last = 16),
#                     collapse = "_"), "_",
#              substring(featureName[1][1], first = 18),
#              "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_v010720.pdf"),
#       plot = bp,
#       height = 8, width = 14)
#} else {
#ggsave(paste0(sub("\\w+\\/$", "", outDir[1]),
#              orderingFactor, "_bargraph_allpops_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePlot),
#              "_in_quantile", quantileNo, "_of_", quantiles,
#              "_by_log2_", libName, "_control_in_", region, "_of_",
#              substring(featureName[1][1], first = 1, last = 5), "_in_",
#              paste0(substring(featureName, first = 10, last = 16),
#                     collapse = "_"), "_",
#              substring(featureName[1][1], first = 18),
#              "_ann_with_GO_BP_enrichment_GO:", GO_ID, "_v010720.pdf"),
#       plot = bp,
#       height = 8, width = 14)
#}
