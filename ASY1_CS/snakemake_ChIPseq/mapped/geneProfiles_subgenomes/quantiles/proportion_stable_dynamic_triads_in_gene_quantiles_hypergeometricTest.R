#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 17.04.2020

# Perform hypergeometric tests to determine whether each ASY1, DMC1 or cM/Mb
# gene quantile is over-represented or under-represented for genes
# assigned to stable, middle, or dynamic homoeolog triad expression bias categories
# (i.e., a stable triad corresponds to a homoeolog triad in the bottom 10%
# in terms of its mean Euclidean distance of the normalised TPM values across 15 tissues
# from the global mean normalised TPM values for that triad;
# and a dynamic triad corresponds to a homoeolog triad in the top 10%
# in terms of its mean Euclidean distance of the normalised TPM values across 15 tissues
# from the global mean normalised TPM values for that triad;
# as defined across 15 tissues in
# Ramirez-Gonazalez et al. 2018 Science 361)
# Also see
# https://github.com/Uauy-Lab/WheatHomoeologExpression/blob/master/02.%20Calculate%20triad%20category.ipynb
# (e.g., is the proportion of genes within a given quantile that
# form part of a dynamic triad
# significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_category) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_category) [m] + ( length(genome_genes) - length(genome_category)) [n]

# Usage 
# ./proportion_stable_dynamic_triads_in_gene_quantiles_hypergeometricTest.R 'ASY1_CS_Rep1_ChIP' 'genes' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome' 100000 6

library(methods)
library(plotrix)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(parallel)

#libName <- "ASY1_CS_Rep1_ChIP"
#featRegion <- "genes"
#quantileFirst <- 1
#quantileLast <- 4
#region <- "genomewide"
#genomeName <- "Agenome_Bgenome_Dgenome"
#samplesNum <- 100000
#minConditions <- 6

args <- commandArgs(trailingOnly = TRUE)
libName <- args[1]
featRegion <- args[2]
quantileFirst <- as.integer(args[3])
quantileLast <- as.integer(args[4])
region <- args[5]
genomeName <- args[6]
samplesNum <- as.numeric(args[7])
minConditions <- as.numeric(args[8])

if(libName %in% "cMMb") {
  outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", featRegion, "/hypergeometricTests/")
}
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
if(libName %in% "cMMb") {
  outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/triad_movement/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", featRegion, "/hypergeometricTests/triad_movement/")
}
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

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
#quantileColours <- makeTransparent(quantileColours)

# Load feature quantiles
if(libName %in% "cMMb") {
  featuresDF <- read.table(paste0(sub("hypergeometricTests/triad_movement/", "", outDir),
                                  "/WesternEurope/features_", quantileLast, "quantiles_by_",
                                  sub("_\\w+$", "", libName), "_of_genes_in_",
                                  genomeName, "_", region, "_WesternEurope.txt"),
                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
} else {
  featuresDF <- read.table(paste0(sub("hypergeometricTests/triad_movement/", "", outDir),
                                  "/WesternEurope/features_", quantileLast, "quantiles_by_",
                                  sub("_\\w+$", "", libName), "_in_", featRegion, "_of_genes_in_",
                                  genomeName, "_", region, "_WesternEurope.txt"),
                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)

}
featuresDF$featureID <- sub(pattern = "\\.\\d+", replacement = "",
                            x = featuresDF$featureID)

# Load table of genes assigned to the different
# homoeolog expression bias categories
triads <- readRDS("/home/ajt200/analysis/wheat/RNAseq_RamirezGonzalez_Uauy_2018_Science/EI_grassroots_data_repo/TablesForExploration/Triads.rds")
triads$gene <- sub(pattern = "1G", replacement = "2G",
                   x = triads$gene)
triads$description <- sub(pattern = "Central", replacement = "Balanced",
                          x = triads$description)
triads$general_description <- sub(pattern = "Central", replacement = "Balanced",
                                  x = triads$general_description)
colnames(triads) <- sub(pattern = "Central", replacement = "Balanced",
                        x = colnames(triads))
triads <- data.frame(triads,
                     description2 = "",
                     stringsAsFactors = F)
triads[(triads$description == "Balanced") , ]$description2 <- "Balanced"
triads[(triads$description == "A.dominant" &
        triads$chr_group == "A") , ]$description2 <- "A.dominant"
triads[(triads$description == "B.dominant" &
        triads$chr_group == "B") , ]$description2 <- "B.dominant"
triads[(triads$description == "D.dominant" &
        triads$chr_group == "D") , ]$description2 <- "D.dominant"
triads[(triads$description == "A.suppressed" &
        triads$chr_group == "A") , ]$description2 <- "A.suppressed"
triads[(triads$description == "B.suppressed" &
        triads$chr_group == "B") , ]$description2 <- "B.suppressed"
triads[(triads$description == "D.suppressed" &
        triads$chr_group == "D") , ]$description2 <- "D.suppressed"
triads[(triads$description == "A.dominant" &
        triads$chr_group != "A") |
       (triads$description == "B.dominant" &
        triads$chr_group != "B") |
       (triads$description == "D.dominant" &
        triads$chr_group != "D") , ]$description2 <- "non_dominant"
triads[(triads$description == "A.suppressed" &
        triads$chr_group != "A") |
       (triads$description == "B.suppressed" &
        triads$chr_group != "B") |
       (triads$description == "D.suppressed" &
        triads$chr_group != "D") , ]$description2 <- "non_suppressed"

## Example using rdist package for calculating for a given triad
## the mean distance across all tissues from the global mean normalised TPM values
## group_id 378 corresponds to the first row (triad) in tm
#group_id <- 378
#local_triad_all_mean <- triads[triads$dataset == "HC_CS_no_stress" &
#                               triads$factor == "all_mean_filter" &
#                               triads$group_id == group_id,]
#local_centroid <- t(matrix(data = local_triad_all_mean$normalised_triad,
#                           nrow = 3,
#                           dimnames = list(local_triad_all_mean$chr_group,
#                                           as.character(unique(local_triad_all_mean$factor)))))
#local_triad_each_factor <- triads[triads$dataset == "HC_CS_no_stress" &
#                                  triads$factor != "all_mean_filter" &
#                                  triads$group_id == group_id,]
#local_mat <- t(matrix(data = local_triad_each_factor$normalised_triad,
#                      nrow = 3,
#                      dimnames = list(local_triad_each_factor$chr_group[1:3],
#                                      as.character(unique(local_triad_each_factor$factor)))))
#library(rdist)
#dists <- cdist(local_mat, local_centroid)
#mean_dist <- mean(dists)

# Subset genes in triads dataframe to only those within a given subgenome
genomeLetter <- unlist(strsplit(gsub("genome", "", genomeName), split = "_"))
if(length(genomeLetter) == 1) {
  triads <- triads[triads$chr_group == genomeLetter,]
}

category_factor <- c("stable", "middle", "dynamic")
category_factor_plot <- c("Stable", "Middle", "Dynamic")
category_colour <- c("Stable 10%" = "dodgerblue3", "Midddle 80%" = "grey60", "Dynamic 10%" = "firebrick1")

# Load table of gene homoeolog triads containing distances
# across tissues from global mean normalised TPM values
triads_movement <- readRDS("/home/ajt200/analysis/wheat/RNAseq_RamirezGonzalez_Uauy_2018_Science/EI_grassroots_data_repo/TablesForExploration/TriadMovement.rds")

# Obtain the names of datasets derived from different tissue types and conditions  
dataset <- unique(triads$dataset)
# Remove "HC_root" (root samples; index 9) and "HC_abiotic_stress_control" (controls for the samples with abiotic stress; index 14)
# datasets as genes in triads are expressed across too few samples (< 6) in these categories
dataset <- dataset[-c(9, 14)]
## Remove "HC_root" (root samples; index 9)
## datasets as genes in triads are expressed across too few samples (< 1) in these categories
#dataset <- dataset[-c(9)]

mclapply(seq_along(dataset), function(xx) {
#for(xx in seq_along(dataset)) {
  print(dataset[xx])
  if(libName %in% "cMMb") {
    outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/triad_movement/", dataset[xx], "/")
  } else {
    outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                     "_in_", featRegion, "/hypergeometricTests/triad_movement/", dataset[xx], "/")
  }
  system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
  plotDir <- paste0(outDir, "plots/")
  system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

  tm <- triads_movement[triads_movement$dataset == dataset[xx] &
                        triads_movement$factor_count >= minConditions,]
  tm_stable <- tm[tm$central_mean_distance <= quantile(tm$central_mean_distance, 0.1),]
  tm_middle <- tm[tm$central_mean_distance > quantile(tm$central_mean_distance, 0.1) &
                  tm$central_mean_distance < quantile(tm$central_mean_distance, 0.9),]
  tm_dynamic <- tm[tm$central_mean_distance >= quantile(tm$central_mean_distance, 0.9),]

  triads_dataset <- triads[triads$dataset == dataset[xx] &
                           triads$factor == "all_mean_filter",]
  triads_dataset <- data.frame(triads_dataset,
                               movement = "",
                               stringsAsFactors = F)
  triads_dataset[triads_dataset$group_id %in% tm_stable$group_id,]$movement <- "stable"
  triads_dataset[triads_dataset$group_id %in% tm_middle$group_id,]$movement <- "middle"
  triads_dataset[triads_dataset$group_id %in% tm_dynamic$group_id,]$movement <- "dynamic"
  triads_dataset <- triads_dataset[triads_dataset$movement != "",]

  genome_genes <- featuresDF$featureID[featuresDF$featureID %in% triads_dataset$gene]
  quantile_genes_list <- lapply(quantileFirst:quantileLast, function(x) {
    featuresDF[featuresDF$quantile == paste0("Quantile ", x) &
               featuresDF$featureID %in% genome_genes,]$featureID
  })
  for(yy in seq_along(category_factor)) {
    print(category_factor[yy])
    genome_category <- as.character(triads_dataset[grepl(pattern = category_factor[yy],
                                                         x = triads_dataset$movement,
                                                         fixed = T),]$gene)
    
    # Get the intersection of genome_category and genome_genes
    # (this excludes genome_category genes not assigned to a chromosome)
    genome_category <- intersect(genome_category, genome_genes)
    
    # Set class for hypergeometric test results object
    setClass("hypergeomTest",
             representation(alternative = "character",
                            alpha0.05 = "numeric",
                            pval = "numeric",
                            observed = "numeric",
                            expected = "numeric",
                            log2obsexp = "numeric",
                            log2alpha = "numeric",
                            quantile_genes = "numeric",
                            proportion_of_quantile = "numeric",
                            random_proportions_of_quantile = "numeric",
                            hypergeomDist = "numeric"))
    
    # P-value is the probability of drawing >= length(quantile_category) [x] features
    # in a sample size of length(quantile_genes) [k] from a total feature set consisting of
    # length(genome_category) [m] + ( length(genome_genes) - length(genome_category)) [n]
    
    # From Karl Broman's answer at
    # https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
    # dhyper(x, m, n, k) gives the probability of drawing exactly x.
    # So P-value is given by the sum of the probabilities of drawing
    # length(quantile_category) to length(quantile_genes)
    
    #lapply(seq_along(quantile_genes_list), function(z) {
    for(z in seq_along(quantile_genes_list)) {
      print(paste0("Quantile ", z))
      quantile_genes <- quantile_genes_list[[z]]
      # Get intersection of gene IDs in quantile z and gene IDs of genome_category
      quantile_category <- intersect(quantile_genes, genome_category)
    
      # Calculate the P-values for over-representation and under-representation
      # of genome_category among quantile z genes
      set.seed(2847502)
      # Over-representation:
      Pval_overrep <- sum(dhyper(x = length(quantile_category):length(quantile_genes),
                                 m = length(genome_category),
                                 n = length(genome_genes) - length(genome_category),
                                 k = length(quantile_genes)))
      print(Pval_overrep)
    
      # Or by 1 minus the sum of the probabilities of drawing 0:(length(quantile_category)-1)
      print(1 - sum(dhyper(x = 0:(length(quantile_category)-1),
                           m = length(genome_category),
                           n = length(genome_genes) - length(genome_category),
                           k = length(quantile_genes))))
    
      # Under-representation
      Pval_underrep <- phyper(q = length(quantile_category),
                              m = length(genome_category),
                              n = length(genome_genes) - length(genome_category),
                              k = length(quantile_genes))
      print(Pval_underrep)
    
      # Sample without replacement
      hgDist <- rhyper(nn = samplesNum,
                       m = length(genome_category),
                       n = length(genome_genes) - length(genome_category),
                       k = length(quantile_genes))
    
      # Calculate P-values and significance levels
      if(length(quantile_category) > mean(hgDist)) {
        Pval <- Pval_overrep
        MoreOrLessThanRandom <- "MoreThanRandom"
        alpha0.05 <- quantile(hgDist, probs = 0.95)[[1]]
      } else {
        Pval <- Pval_underrep
        MoreOrLessThanRandom <- "LessThanRandom"
        alpha0.05 <- quantile(hgDist, probs = 0.05)[[1]]
      }
    
      hgTestResults <- new("hypergeomTest",
                           alternative = MoreOrLessThanRandom,
                           alpha0.05 = alpha0.05,
                           pval = Pval,
                           observed = length(quantile_category),
                           expected = mean(hgDist),
                           log2obsexp = log2( length(quantile_category) / mean(hgDist) ),
                           log2alpha  = log2( alpha0.05 / mean(hgDist) ),
                           quantile_genes = length(quantile_genes),
                           proportion_of_quantile = length(quantile_category) / length(quantile_genes),
                           random_proportions_of_quantile = hgDist / length(quantile_genes),
                           hypergeomDist = hgDist)
      if(libName %in% "cMMb") {
      save(hgTestResults,
           file = paste0(outDir,
                         dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                         "_by_", libName, "_of_genes_in_",
                         genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      } else {
      save(hgTestResults,
           file = paste0(outDir,
                         dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                         "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                         genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      }
    
      # Generate histogram
      if(libName %in% "cMMb") {
      pdf(paste0(plotDir,
                 dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                 "_by_", libName, "_of_genes_in_",
                 genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, "_hist.pdf"),
                 height = 4.5, width = 5)
      } else {
      pdf(paste0(plotDir,
                 dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                 "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                 genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, "_hist.pdf"),
                 height = 4.5, width = 5)
      }
      par(mar = c(3.1, 3.1, 4.1, 1.1),
          mgp = c(1.85, 0.75, 0))
      ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
      #options(scipen = 100)
      # Calculate max density
      maxDensityPlus <- max(density(hgTestResults@hypergeomDist)$y)*1.2
      if(hgTestResults@alternative == "MoreThanRandom") {
        xlim <- c(pmin(min(hgTestResults@hypergeomDist)/1.2),
                  pmax(hgTestResults@observed*1.2, hgTestResults@alpha0.05*1.2))
        textX1 <- quantile(xlim, 0.25)[[1]]
    #    textX1 <- min(hgTestResults@hypergeomDist)/1.15
      } else {
        xlim <- c(pmin(hgTestResults@observed/1.2),
                  max(hgTestResults@hypergeomDist)*1.2)
        textX1 <- quantile(xlim, 0.75)[[1]]
    #    textX1 <- min(hgTestResults@hypergeomDist)/1.15
      }
      hist(hgTestResults@hypergeomDist,
           breaks = 50,
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
           at = pretty(density(hgTestResults@hypergeomDist)$y),
           lwd = 2)
      mtext(side = 2,
            text = "Density",
            line = 1.85)
      axis(side = 1,
           at = pretty(xlim),
           lwd = 2)
      mtext(side = 1,
            text = bquote("Genes"),
            line = 1.85)
      titleText <- list(bquote(.(category_factor[yy]) ~ "genes in" ~
                               .(sub("_\\w+$", "", libName)) ~ "Quantile" ~ .(as.character(z)) ~
                               "(" * .(featRegion) * ") in" ~
                               .(gsub("_", " ", genomeName)) ~ .(region)),
                        bquote(italic("P")*" = "*
    #                           .(as.character(round(hgTestResults@pval,
    #                                                digits = 6)))),
                               .(as.character(hgTestResults@pval))),
                        bquote("Samples (hypergeometric distribution) = "*.(prettyNum(length(hgTestResults@hypergeomDist),
                                                                                      big.mark = ",",
                                                                                      trim = T))))
      mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = c(0.7, 1, 1))
      lines(density(hgTestResults@hypergeomDist),
            col = "dodgerblue3",
            lwd = 1.5)
      ablineclip(v = hgTestResults@expected,
                 y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
      ablineclip(v = hgTestResults@observed,
                 y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
      ablineclip(v = hgTestResults@alpha0.05,
                 y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
      text(x = c(textX1,
                 hgTestResults@expected,
                 hgTestResults@observed,
                 hgTestResults@alpha0.05),
           y = c(maxDensityPlus*.95,
                 maxDensityPlus,
                 maxDensityPlus,
                 maxDensityPlus*.95),
           labels = c("Simulated",
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
    
    
    options(scipen = 100)
    
    # Plot bar graph summarising permutation test results
    pt_list <- list()
    for(z in quantileFirst:quantileLast) {
      if(libName %in% "cMMb") {
      load(paste0(outDir,
                  dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_", libName, "_of_genes_in_",
                  genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      } else {
      load(paste0(outDir,
                  dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                  genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
     }
      pt_list <- c(pt_list, hgTestResults)
    }
    bargraph_df <- data.frame(Quantile = paste0("Quantile ", quantileFirst:quantileLast),
                              log2ObsExp = sapply(seq_along(pt_list), function(x) { pt_list[[x]]@log2obsexp }),
                              log2alpha0.05 = sapply(seq_along(pt_list), function(x) { pt_list[[x]]@log2alpha }))
    bargraph_df$Quantile <- factor(bargraph_df$Quantile,
                                   levels = paste0("Quantile ", quantileFirst:quantileLast))
    bp <- ggplot(data = bargraph_df,
                 mapping = aes(x = Quantile,
                               y = log2ObsExp,
                               fill = " ")) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "",
                        values = c("dodgerblue3"),
                        labels = " ") +
      geom_point(mapping = aes(x = Quantile,
                               y = log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "grey80", size = 20) +
      labs(y = bquote("Log"[2]*"(observed/expected) genes in quantile")) +
    #  scale_y_continuous(limits = c(-1.5, 1.5)) +
      scale_x_discrete(position = "top") +
      guides(fill = guide_legend(direction = "horizontal",
                                 label.position = "top",
                                 label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
                                 nrow = 1,
                                 byrow = TRUE)) +
      theme_bw() +
      theme(axis.line.y = element_line(size = 1, colour = "black"),
            axis.ticks.y = element_line(size = 1, colour = "black"),
            axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
            axis.title.y = element_text(size = 20, colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, angle = 90),
            axis.title.x = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            #legend.position = c(0.05, 0.30),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(colour = "transparent",
                                      fill = "transparent"),
            plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
            plot.title = element_text(size = 16, colour = "black", hjust = 0.5)) +
      ggtitle(bquote(.(category_factor[yy]) ~ "genes in" ~
                     .(sub("_\\w+$", "", libName)) ~ "quantiles" ~
                     "(" * .(featRegion) * ") in" ~
                     .(gsub("_", " ", genomeName)) ~ .(region) ~
                     "(" * .(prettyNum(samplesNum,
                                       big.mark = ",",
                                       trim = T)) ~ " samples)"))
    if(libName %in% "cMMb") {
    ggsave(paste0(plotDir,
                  "bargraph_", dataset[xx], "_", category_factor[yy], "_representation_among_", quantileLast,
                  "quantiles_by_", libName, "_of_genes_in_",
                  genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, ".pdf"),
           plot = bp,
           height = 8, width = 12)
    } else {
    ggsave(paste0(plotDir,
                  "bargraph_", dataset[xx], "_", category_factor[yy], "_representation_among_", quantileLast,
                  "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                  genomeName, "_", region, "_hypergeomTestRes_minConditions", minConditions, ".pdf"),
           plot = bp,
           height = 8, width = 12)
    }
  }
#}
}, mc.cores = length(dataset), mc.preschedule = F)
