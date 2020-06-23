#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 15.04.2020

# For all three wheat subgenomes, load and plot as bargraphs results from
# hypergeometric tests to determine whether each ASY1, DMC1 or cM/Mb
# gene quantile is over-represented or under-represented for
# genes assigned to homoeolog expression bias categories
# (i.e., "Balanced",# ".dominant", "non_dominant", ".suppressed", "non_suppressed",
# "A.dominant", "B.dominant", "D.dominant",
# "A.suppressed", "B.suppressed", "D.suppressed";
# as defined for various tissue types and conditions in
# Ramirez-Gonazalez et al. 2018 Science 361)
# (e.g., is the proportion of genes within a given quantile that
# exhibit "A.dominant" expression bias
# significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_category) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_category) [m] + ( length(genome_genes) - length(genome_category)) [n]

# Usage 
# ./proportion_homoeolog_expression_bias_triad_movement_categories_in_gene_quantiles_hypergeometricTest_bargraph_only.R 'ASY1_CS_Rep1_ChIP' 'genes' 1 4 'genomewide' Agenome_Bgenome_Dgenome 100000 6

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
  outDir1 <- paste0("quantiles_by_", libName, "/hypergeometricTests/homoeolog_exp_bias/")
  outDir2 <- paste0("quantiles_by_", libName, "/hypergeometricTests/triad_movement/")
} else {
  outDir1 <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                    "_in_", featRegion, "/hypergeometricTests/homoeolog_exp_bias/")
  outDir2 <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                    "_in_", featRegion, "/hypergeometricTests/triad_movement/")
}
dataset1 <- system(paste0("ls ", outDir1), intern = T)
dataset2 <- system(paste0("ls ", outDir2), intern = T)
dataset <- intersect(dataset1, dataset2)

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

# Subset genes in triads dataframe to only those within a given subgenome
genomeLetter <- unlist(strsplit(gsub("genome", "", genomeName), split = "_"))
if(length(genomeLetter) == 1) {
  category_factor1 <- c("Balanced",
                       ".dominant", "non_dominant", ".suppressed", "non_suppressed")
  category_factor_plot1 <- c("Balanced",
                            "Dominant", "Non-dominant", "Suppressed", "Non-suppressed")
  category_colour1 <- c("Balanced" = "#AAAAAA",
                       "Dominant" = "turquoise4", "Non-dominant" = "paleturquoise3", "Suppressed" = "lightgoldenrod", "Non-suppressed" = "goldenrod4")
} else {
  category_factor1 <- c("Balanced",
                       ".dominant", "non_dominant", ".suppressed", "non_suppressed",
                       "A.dominant", "B.dominant", "D.dominant",
                       "A.suppressed", "B.suppressed", "D.suppressed")
  category_factor_plot1 <- c("Balanced",
                            "Dominant", "Non-dominant", "Suppressed", "Non-suppressed",
                            "A dominant", "B dominant", "D dominant",
                            "A suppressed", "B suppressed", "D suppressed")
  category_colour1 <- c("Balanced" = "#AAAAAA",
                       "Dominant" = "turquoise4", "Non-dominant" = "paleturquoise3", "Suppressed" = "lightgoldenrod", "Non-suppressed" = "goldenrod4",
                       "A dominant" = "#579D1C", "B dominant" = "#4B1F6F", "D dominant" = "#FF950E",
                       "A suppressed" = "#b2E08a", "B suppressed" = "#0eC7ff", "D suppressed" ="#ffCF0e")
}

category_factor2 <- c("stable", "middle", "dynamic")
category_factor_plot2 <- c("Stable 10%", "Middle 80%", "Dynamic 10%")
category_colour2 <- c("Stable 10%" = "green3", "Middle 80%" = "grey60", "Dynamic 10%" = "deeppink3")

category_factor <- c(category_factor1, category_factor2)
category_factor_plot <- c(category_factor_plot1, category_factor_plot2)
category_colour <- c(category_colour1, category_colour2)

options(scipen = 100)

mclapply(seq_along(dataset), function(xx) {
#for(xx in seq_along(dataset)) {
  if(libName %in% "cMMb") {
    outDir1 <- paste0("quantiles_by_", libName, "/hypergeometricTests/homoeolog_exp_bias/", dataset[xx], "/")
    outDir2 <- paste0("quantiles_by_", libName, "/hypergeometricTests/triad_movement/", dataset[xx], "/")
  } else {
    outDir1 <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                      "_in_", featRegion, "/hypergeometricTests/homoeolog_exp_bias/", dataset[xx], "/")
    outDir2 <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                      "_in_", featRegion, "/hypergeometricTests/triad_movement/", dataset[xx], "/")
  }
  plotDir <- paste0(outDir1, "plots/combined_bargraph/")
  system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
  
  # Plot bar graph summarising permutation test results,
  # including subgenome-specific expression bias categories
  hg_list1 <- lapply(seq_along(category_factor1), function(y) {
    hg_list_quantile <- list() 
    for(z in quantileFirst:quantileLast) {
      if(libName %in% "cMMb") {
      load(paste0(outDir1,
                  dataset[xx], "_", category_factor1[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_", libName, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes.RData"))
      } else {
      load(paste0(outDir1,
                  dataset[xx], "_", category_factor1[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes.RData"))
      }
      hg_list_quantile <- c(hg_list_quantile, hgTestResults)
    }
    return(hg_list_quantile)
  })
  hg_list2 <- lapply(seq_along(category_factor2), function(y) {
    hg_list_quantile <- list() 
    for(z in quantileFirst:quantileLast) {
      if(libName %in% "cMMb") {
      load(paste0(outDir2,
                  dataset[xx], "_", category_factor2[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_", libName, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      } else {
      load(paste0(outDir2,
                  dataset[xx], "_", category_factor2[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      }
      hg_list_quantile <- c(hg_list_quantile, hgTestResults)
    }
    return(hg_list_quantile)
  })
  hg_list <- c(hg_list1, hg_list2)
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
                               label.theme = element_text(size = 30, hjust = 0, vjust = 0.5, angle = 90),
                               nrow = 3,
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
          axis.text.y = element_text(size = 40, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
          axis.title.y = element_text(size = 40, colour = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 40, colour = "black", hjust = 0.5, vjust = 0.5, angle = 0),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
#          legend.position = "none",
#          legend.position = c(0.05, 0.30),
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(colour = "transparent",
                                    fill = "transparent"),
          plot.margin = unit(c(5.5, 5.5, 40.5, 5.5), "pt"),
          plot.title = element_text(size = 40, colour = "black", hjust = 0.5)) +
    ggtitle(bquote(.(dataset[xx]) ~ "homoeolog expression categories in" ~
                   .(sub("_\\w+$", "", libName)) ~ "quantiles" ~
                   "(" * .(featRegion) * ")" ~
                   "(" * .(prettyNum(samplesNum,
                                     big.mark = ",",
                                     trim = T)) ~ "samples)"))
  if(libName %in% "cMMb") {
  ggsave(paste0(plotDir,
                "combined_bargraph_", dataset[xx], "_homoeolog_expresssion_variation_category_representation_among_", quantileLast,
                "quantiles_by_", libName, "_of_genes_in_", genomeName, "_",
                region, "_hypergeomTestRes.pdf"),
         plot = bp,
         height = 12.5, width = 50, limitsize = F)
  } else {
  ggsave(paste0(plotDir,
                "combined_bargraph_", dataset[xx], "_homoeolog_expresssion_variation_category_representation_among_", quantileLast,
                "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                region, "_hypergeomTestRes.pdf"),
         plot = bp,
         height = 12.5, width = 50, limitsize = F)
  }
  # Plot bar graph summarising permutation test results,
  # excluding subgenome-specific expression bias categories
  if(length(genomeLetter) > 1) {
  category_factor <- category_factor[-c(6:11)]
  category_factor_plot <- category_factor_plot[-c(6:11)]
  category_colour <- category_colour[-c(6:11)]
  category_factor1 <- category_factor1[-c(6:11)]
  category_factor_plot1 <- category_factor_plot1[-c(6:11)]
  category_colour1 <- category_colour1[-c(6:11)]
  hg_list1 <- lapply(seq_along(category_factor1), function(y) {
    hg_list_quantile <- list() 
    for(z in quantileFirst:quantileLast) {
      if(libName %in% "cMMb") {
      load(paste0(outDir1,
                  dataset[xx], "_", category_factor1[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_", libName, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes.RData"))
      } else {
      load(paste0(outDir1,
                  dataset[xx], "_", category_factor1[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes.RData"))
      }
      hg_list_quantile <- c(hg_list_quantile, hgTestResults)
    }
    return(hg_list_quantile)
  })
  hg_list2 <- lapply(seq_along(category_factor2), function(y) {
    hg_list_quantile <- list() 
    for(z in quantileFirst:quantileLast) {
      if(libName %in% "cMMb") {
      load(paste0(outDir2,
                  dataset[xx], "_", category_factor2[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_", libName, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      } else {
      load(paste0(outDir2,
                  dataset[xx], "_", category_factor2[y], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                  region, "_hypergeomTestRes_minConditions", minConditions, ".RData"))
      }
      hg_list_quantile <- c(hg_list_quantile, hgTestResults)
    }
    return(hg_list_quantile)
  })
  hg_list <- c(hg_list1, hg_list2)
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
          axis.text.x = element_text(size = 32, colour = "black", hjust = 0.5, vjust = 0.5, angle = 0),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
#          legend.position = "none",
#          legend.position = c(0.05, 0.30),
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(colour = "transparent",
                                    fill = "transparent"),
          plot.margin = unit(c(5.5, 5.5, 40.5, 5.5), "pt"),
          plot.title = element_text(size = 24, colour = "black", hjust = 0.5)) +
    ggtitle(bquote(.(dataset[xx]) ~ "homoeolog expression categories in" ~
                   .(sub("_\\w+$", "", libName)) ~ "quantiles" ~
                   "(" * .(featRegion) * ")" ~
                   "(" * .(prettyNum(samplesNum,
                                     big.mark = ",",
                                     trim = T)) ~ "samples)"))
  if(libName %in% "cMMb") {
  ggsave(paste0(plotDir,
                "combined_bargraph_", dataset[xx], "_homoeolog_expresssion_variation_category_representation_among_", quantileLast,
                "quantiles_by_", libName, "_of_genes_in_", genomeName, "_",
                region, "_hypergeomTestRes_excl_subgenomes.pdf"),
         plot = bp,
         height = 8, width = 32, limitsize = F)
  } else {
  ggsave(paste0(plotDir,
                "combined_bargraph_", dataset[xx], "_homoeolog_expresssion_variation_category_representation_among_", quantileLast,
                "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_", genomeName, "_",
                region, "_hypergeomTestRes_excl_subgenomes.pdf"),
         plot = bp,
         height = 8, width = 32, limitsize = F)
  }
  }
#}
}, mc.cores = length(dataset), mc.preschedule = F)
