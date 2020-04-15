#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 15.04.2020

# Perform hypergeometric tests to determine whether each ASY1, DMC1 or cM/Mb
# gene quantile is over-represented or under-represented for
# genes assigned to homoeolog expression bias categories
# (i.e., "Balanced", "A dominant", "B dominant", "D dominant",
# "A suppressed", "B suppressed", "D suppressed";
# as defined for various tissue types and conditions in
# Ramirez-Gonazalez et al. 2018 Science 361) 
# (e.g., is the proportion of genes within a given quantile that
# are within a triad that exhibits "A dominant" expression bias
# significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_category) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_category) [m] + ( length(genome_genes) - length(genome_category)) [n]

# Usage 
# ./proportion_homoeolog_expression_bias_categories_in_gene_quantiles_hypergeometricTest.R 'ASY1_CS_Rep1_ChIP' 'genes' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome' 100000

library(methods)
library(plotrix)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

#libName <- "ASY1_CS_Rep1_ChIP"
#featRegion <- "genes"
#quantileFirst <- 1
#quantileLast <- 4
#region <- "genomewide"
#genomeName <- "Agenome_Bgenome_Dgenome"
#samplesNum <- 100000

args <- commandArgs(trailingOnly = TRUE)
libName <- args[1]
featRegion <- args[2]
quantileFirst <- as.integer(args[3])
quantileLast <- as.integer(args[4])
region <- args[5]
genomeName <- args[6]
samplesNum <- as.numeric(args[7])

if(libName %in% "cMMb") {
  outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", featRegion, "/hypergeometricTests/")
}
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
if(libName %in% "cMMb") {
  outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/homoeolog_exp_bias/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", featRegion, "/hypergeometricTests/homoeolog_exp_bias/")
}
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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

# Load feature quantiles
if(libName %in% "cMMb") {
  featuresDF <- read.table(paste0(sub("hypergeometricTests/homoeolog_exp_bias/", "", outDir),
                                  "/WesternEurope/features_", quantileLast, "quantiles_by_",
                                  sub("_\\w+$", "", libName), "_of_genes_in_",
                                  genomeName, "_", region, "_WesternEurope.txt"),
                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
} else {
  featuresDF <- read.table(paste0(sub("hypergeometricTests/homoeolog_exp_bias/", "", outDir),
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

# Subset genes in triads dataframe to only those within a given subgenome
genomeLetter <- unlist(strsplit(gsub("genome", "", genomeName), split = "_"))
if(length(genomeLetter) == 1) {
  triads <- triads[triads$chr_group == genomeLetter,]
  category_factor <- c("Balanced",
                       ".dominant", "non_dominant", ".suppressed", "non_suppressed")
  category_colour <- c(Balanced = "#AAAAAA",
                       .dominant = "darkcyan", non_dominant = "cyan", .suppressed = "goldenrod", non_suppressed = "goldenrod4")
} else {
  category_factor <- c("Balanced",
                       ".dominant", "non_dominant", ".suppressed", "non_suppressed",
                       "A.dominant", "B.dominant", "D.dominant",
                       "A.suppressed", "B.suppressed", "D.suppressed")
  category_colour <- c(Balanced = "#AAAAAA",
                       .dominant = "darkcyan", non_dominant = "cyan", .suppressed = "goldenrod", non_suppressed = "goldenrod4",
                       A.dominant = "#579D1C", B.dominant = "#4B1F6F", D.dominant = "#FF950E",
                       A.suppressed = "#b2E08a", B.suppressed = "#0eC7ff", D.suppressed ="#ffCF0e")
}

# Obtain the names of datasets derived from different tissue types and conditions  
dataset <- unique(triads$dataset)

for(xx in seq_along(dataset)) {
  triads_dataset <- triads[triads$dataset == dataset[xx] &
                           triads$factor == "all_mean_filter",]
  genome_genes <- featuresDF$featureID[featuresDF$featureID %in% triads_dataset$gene]
  quantile_genes_list <- lapply(quantileFirst:quantileLast, function(x) {
    featuresDF[featuresDF$quantile == paste0("Quantile ", x) &
               featuresDF$featureID %in% genome_genes,]$featureID
  })
  for(yy in seq_along(category_factor)) {
    genome_category <- as.character(triads_dataset[grepl(pattern = category_factor[yy],
                                                         x = triads_dataset$description2,
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
                         genomeName, "_", region, "_hypergeomTestRes.RData"))
      } else {
      save(hgTestResults,
           file = paste0(outDir,
                         dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                         "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                         genomeName, "_", region, "_hypergeomTestRes.RData"))
      }
    
      # Generate histogram
      if(libName %in% "cMMb") {
      pdf(paste0(plotDir,
                 dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                 "_by_", libName, "_of_genes_in_",
                 genomeName, "_", region, "_hypergeomTestRes_hist.pdf"),
                 height = 4.5, width = 5)
      } else {
      pdf(paste0(plotDir,
                 dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                 "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                 genomeName, "_", region, "_hypergeomTestRes_hist.pdf"),
                 height = 4.5, width = 5)
      }
      par(mar = c(3.1, 3.1, 4.1, 1.1),
          mgp = c(1.85, 0.75, 0))
      ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
      #options(scipen = 100)
      # Calculate max density
      maxDensityPlus <- max(density(hgTestResults@hypergeomDist)$y)*1.2
      if(hgTestResults@alternative == "MoreThanRandom") {
        xlim <- c(pmin(0, min(hgTestResults@hypergeomDist)/1.2),
                  pmax(hgTestResults@observed*1.2, hgTestResults@alpha0.05*1.2))
        textX1 <- quantile(xlim, 0.25)[[1]]
    #    textX1 <- min(hgTestResults@hypergeomDist)/1.15
      } else {
        xlim <- c(pmin(0, hgTestResults@observed/1.2),
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
                  genomeName, "_", region, "_hypergeomTestRes.RData"))
      } else {
      load(paste0(outDir,
                  dataset[xx], "_", category_factor[yy], "_representation_among_quantile", z, "_of_", quantileLast,
                  "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                  genomeName, "_", region, "_hypergeomTestRes.RData"))
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
                  genomeName, "_", region, "_hypergeomTestRes.pdf"),
           plot = bp,
           height = 8, width = 12)
    } else {
    ggsave(paste0(plotDir,
                  "bargraph_", dataset[xx], "_", category_factor[yy], "_representation_among_", quantileLast,
                  "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                  genomeName, "_", region, "_hypergeomTestRes.pdf"),
           plot = bp,
           height = 8, width = 12)
    }
  }
}
