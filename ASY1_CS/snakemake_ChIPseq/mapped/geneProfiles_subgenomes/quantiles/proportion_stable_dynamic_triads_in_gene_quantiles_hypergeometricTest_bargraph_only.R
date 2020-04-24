#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 17.04.2020

# For all three wheat subgenomes, load and plot as bargraphs results from
# hypergeometric tests to determine whether each ASY1, DMC1 or cM/Mb
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
# ./proportion_stable_dynamic_triads_in_gene_quantiles_hypergeometricTest_bargraph_only.R 'ASY1_CS_Rep1_ChIP' 'genes' 1 4 'genomewide' Agenome_Bgenome_Dgenome 100000 6

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
if(libName %in% "cMMb") {
  outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/triad_movement/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", featRegion, "/hypergeometricTests/triad_movement/")
}
dataset <- system(paste0("ls ", outDir), intern = T)
## Remove "HC_abiotic_stress_control" (controls for the samples with abiotic stress; index 5)
## datasets as genes in triads are expressed across too few samples (< 6) in these categories
## NOTE: this applies only to ASY1 and DMC1 in promoters where I was interested in
## evaluating triad category representation among gene quantiles while imposing a lower
## minimum number of conditions in which a gene triad must be expressed 
#dataset <- dataset[-c(5)]

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

category_factor <- c("stable", "middle", "dynamic")
category_factor_plot <- c("Stable 10%", "Middle 80%", "Dynamic 10%")
category_colour <- c("Stable 10%" = "green3", "Middle 80%" = "grey60", "Dynamic 10%" = "deeppink3")

options(scipen = 100)

mclapply(seq_along(dataset), function(xx) {
#for(xx in seq_along(dataset)) {
  if(libName %in% "cMMb") {
    outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/triad_movement/", dataset[xx], "/")
  } else {
    outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                     "_in_", featRegion, "/hypergeometricTests/triad_movement/", dataset[xx], "/")
  }
  plotDir <- paste0(outDir, "plots/combined_bargraph/")
  system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
  
  # Plot bar graph summarising permutation test results
  hg_list <- lapply(seq_along(category_factor), function(y) {
    hg_list_quantile <- list() 
    for(z in quantileFirst:quantileLast) {
      if(libName %in% "cMMb") {
      load(paste0(outDir,
                  dataset[xx], "_", category_factor[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_", libName, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      } else {
      load(paste0(outDir,
                  dataset[xx], "_", category_factor[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      }
      hg_list_quantile <- c(hg_list_quantile, hgTestResults)
    }
    return(hg_list_quantile)
  })
  bargraph_df <- data.frame(Category = rep(category_factor_plot, each = quantileLast),
                            Quantile = rep(paste0("Quantile ", quantileFirst:quantileLast), times = length(category_factor_plot)),
                            log2ObsExp = c(sapply(seq_along(category_factor), function(y) {
                                             sapply(seq_along(hg_list[[y]]), function(z) {
                                               hg_list[[y]][[z]]@log2obsexp
                                             })
                                           })),
                            log2alpha0.05 = c(sapply(seq_along(category_factor), function(y) {
                                                sapply(seq_along(hg_list[[y]]), function(z) {
                                                  hg_list[[y]][[z]]@log2alpha
                                                })
                                              })))
  bargraph_df$Quantile <- factor(bargraph_df$Quantile,
                                 levels = paste0("Quantile ", quantileFirst:quantileLast))
  bargraph_df$Category <- factor(bargraph_df$Category,
                                 levels = category_factor_plot)
          
  bp <- ggplot(data = bargraph_df,
               mapping = aes(x = Quantile,
                             y = log2ObsExp,
                             fill = Category)) +
    geom_bar(stat = "identity",
             position = position_dodge()) +
    scale_fill_manual(name = "",
                      values = category_colour,
                      labels = levels(bargraph_df$Category)) +
    geom_point(mapping = aes(x = Quantile,
                             y = log2alpha0.05),
               position = position_dodge(0.9),
               shape = "-", colour  = "grey40", size = 20) +
    guides(fill = guide_legend(direction = "horizontal",
                               label.position = "top",
                               label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
                               nrow = 1,
                               byrow = TRUE)) +
    geom_segment(mapping = aes(x = 0.55, y = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1,
                               xend = 1.45, yend = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1),
                 colour = quantileColours[1],
                 inherit.aes = F, size = 5) +
    geom_segment(mapping = aes(x = 1.55, y = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1,
                               xend = 2.45, yend = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1),
                 colour = quantileColours[2],
                 inherit.aes = F, size = 5) +
    geom_segment(mapping = aes(x = 2.55, y = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1,
                               xend = 3.45, yend = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1),
                 colour = quantileColours[3],
                 inherit.aes = F, size = 5) +
    geom_segment(mapping = aes(x = 3.55, y = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1,
                               xend = 4.45, yend = min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1),
                 colour = quantileColours[4],
                 inherit.aes = F, size = 5) +
#    geom_boxplot(data = data.frame(quantile_lines = factor(levels(bargraph_df$Quantile),
#                                                                  levels = levels(bargraph_df$Quantile)),
#                                   yval = rep(min(c(bargraph_df$log2ObsExp, bargraph_df$log2alpha0.05))-0.1, times = quantileLast)),
#                 mapping = aes(x = quantile_lines,
#                               y = yval,
#                               colour = quantile_lines),
#                 inherit.aes = F, size = 3) +
#    scale_colour_manual(name = "",
#                        values = quantileColours,
#                        labels = "") +
    labs(y = bquote("Log"[2]*"(observed/expected) genes in quantile")) +
#    scale_y_continuous(limits = c(-0.8, 0.8)) +
    scale_x_discrete(position = "bottom") +
    theme_bw() +
    theme(axis.line.y = element_line(size = 1, colour = "black"),
          axis.ticks.y = element_line(size = 1, colour = "black"),
          axis.text.y = element_text(size = 25, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
          axis.title.y = element_text(size = 25, colour = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 28, colour = "black", hjust = 0.5, vjust = 0.5, angle = 180),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #legend.position = "none",
          #legend.position = c(0.05, 0.30),
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(colour = "transparent",
                                    fill = "transparent"),
          plot.margin = unit(c(5.5, 5.5, 40.5, 5.5), "pt"),
          plot.title = element_text(size = 14, colour = "black", hjust = 0.5)) +
    ggtitle(bquote(.(dataset[xx]) ~ "triad movement categories in" ~
                   .(sub("_\\w+$", "", libName)) ~ "quantiles" ~
                   "(" * .(featRegion) * ")" ~
                   "(" * .(prettyNum(samplesNum,
                                     big.mark = ",",
                                     trim = T)) ~ "samples)"))
  if(libName %in% "cMMb") {
  ggsave(paste0(plotDir,
                "combined_bargraph_", dataset[xx], "_stable_dynamic_triad_representation_among_", quantileLast,
                "quantiles_by_", libName, "_of_genes_in_", genomeName, "_",
                region, "_hypergeomTestRes_minConditions", minConditions, ".pdf"),
         plot = bp,
         height = 8, width = 12)
  } else {
  ggsave(paste0(plotDir,
                "combined_bargraph_", dataset[xx], "_stable_dynamic_triad_representation_among_", quantileLast,
                "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                region, "_hypergeomTestRes_minConditions", minConditions, ".pdf"),
         plot = bp,
         height = 8, width = 12)
  }
#}
}, mc.cores = length(dataset), mc.preschedule = F)
