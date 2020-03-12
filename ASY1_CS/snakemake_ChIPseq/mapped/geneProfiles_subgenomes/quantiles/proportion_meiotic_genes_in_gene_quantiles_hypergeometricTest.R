#!/applications/R/R-3.5.0/bin/Rscript

# Perform hypergeometric tests to determine whether each ASY1 or DMC1
# gene quantile is over-represented or under-represented for
# Meiotic genes (as defined by Ksenia Krasileva) and, separately,
# meiotic-module genes (as defined in Alabdullah et al. 2019 Front. Plant Sci.)
# (i.e., is the propotion of Meiotic genes or meiotic-module genes
# contained within each gene quantile significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_meioOrthologs) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_meioOrthologs) [m] + ( length(genome_genes) - length(genome_meioOrthologs)) [n]

### Perform hypergeometric test to determine whether a significant proportion of 
### elements belonging to each transposon family are differentially expressed in kyp suvh5 suvh6

# P-value is the probability of drawing >= length(DE_TEfamIDs) [x] TEs
# in a sample size of length(DE_TEIDs) [k] from a total TE set consisting of
# length(TAIR10_TEfamIDs) [m] + (length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs)) [n]

# Usage 
# ./proportion_meiotic_genes_in_gene_quantiles_hypergeometricTest.R 'ASY1_CS_Rep1_ChIP' 'bodies' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome' 100000

library(methods)
library(plotrix)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)


#libName <- "ASY1_CS_Rep1_ChIP"
#featRegion <- "bodies"
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

outDir <- paste0("quantiles_by_", sub("_\\w+$", "", libName), "_in_", featRegion, "/hypergeometricTests/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
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
featuresDF <- read.table(paste0("quantiles_by_", sub("_\\w+$", "", libName), "_in_", featRegion,
                                "/WesternEurope/features_", quantileLast, "quantiles_by_",
                                sub("_\\w+$", "", libName), "_in_", featRegion, "_of_genes_in_",
                                genomeName, "_", region, "_WesternEurope.txt"),
                         header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
featuresDF$featureID <- sub(pattern = "\\.\\d+", replacement = "",
                            x = featuresDF$featureID)
genome_genes <- featuresDF$featureID
quantile_genes_list <- lapply(quantileFirst:quantileLast, function(x) {
  featuresDF[featuresDF$quantile == paste0("Quantile ", x),]$featureID
})
rm(featuresDF); gc()

# Load meio
# Note: these two sets of meiotic genes share 271 common genes that are assigned to a chromosome
meio1 <- read.table("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/Table_S4_meiotic_GO_genes.tsv",
                    header = T, stringsAsFactors = F)
genome_meio1 <- as.character(meio1$Gene.ID)
meio2 <- read.table("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/Table_S4_meiotic_gene_orthologs.tsv",
                    header = T, sep = "\t", stringsAsFactors = F)
genome_meio2 <- as.character(meio2$Gene.ID)

genome_meio <- union(genome_meio1, genome_meio2)

# Get the intersection of genome_meio and genome_genes
# (this excludes meio not assigned to a chromosome)
genome_meio <- intersect(genome_meio, genome_genes)

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

# P-value is the probability of drawing >= length(quantile_meio) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_meio) [m] + ( length(genome_genes) - length(genome_meio)) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(quantile_meio) to length(quantile_genes)

#lapply(seq_along(quantile_genes_list), function(z) {
for(z in seq_along(quantile_genes_list)) {
  quantile_genes <- quantile_genes_list[[z]]
  # Get intersection of gene IDs in quantile z and gene IDs of meio
  quantile_meio <- intersect(quantile_genes, genome_meio)

  # Calculate the P-values for over-representation and under-representation
  # of meio among quantile z genes
  set.seed(2847502)
  # Over-representation:
  Pval_overrep <- sum(dhyper(x = length(quantile_meio):length(quantile_genes),
                             m = length(genome_meio),
                             n = length(genome_genes) - length(genome_meio),
                             k = length(quantile_genes)))
  print(Pval_overrep)

  # Or by 1 minus the sum of the probabilities of drawing 0:(length(quantile_meio)-1)
  print(1 - sum(dhyper(x = 0:(length(quantile_meio)-1),
                       m = length(genome_meio),
                       n = length(genome_genes) - length(genome_meio),
                       k = length(quantile_genes))))

  # Under-representation
  Pval_underrep <- phyper(q = length(quantile_meio),
                          m = length(genome_meio),
                          n = length(genome_genes) - length(genome_meio),
                          k = length(quantile_genes))
  print(Pval_underrep)

  # Sample without replacement
  hgDist <- rhyper(nn = samplesNum,
                   m = length(genome_meio),
                   n = length(genome_genes) - length(genome_meio),
                   k = length(quantile_genes))

  # Calculate P-values and significance levels
  if(length(quantile_meio) > mean(hgDist)) {
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
                       observed = length(quantile_meio),
                       expected = mean(hgDist),
                       log2obsexp = log2( length(quantile_meio) / mean(hgDist) ),
                       log2alpha  = log2( alpha0.05 / mean(hgDist) ),
                       quantile_genes = length(quantile_genes),
                       proportion_of_quantile = length(quantile_meio) / length(quantile_genes),
                       random_proportions_of_quantile = hgDist / length(quantile_genes),
                       hypergeomDist = hgDist)
  save(hgTestResults,
       file = paste0(outDir,
                     "meiotic_gene_representation_among_quantile", z, "_of_", quantileLast,
                     "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
                     genomeName, "_", region, "_hypergeomTestRes.RData"))

  # Generate histogram
  pdf(paste0(plotDir,
             "meiotic_gene_representation_among_quantile", z, "_of_", quantileLast,
             "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
             genomeName, "_", region, "_hypergeomTestRes_hist.pdf"),
             height = 4.5, width = 5)
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
  titleText <- list(bquote("Meiotic genes in" ~
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
  load(
       paste0(outDir,
              "meiotic_gene_representation_among_quantile", z, "_of_", quantileLast,
              "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.RData"))
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
  ggtitle(bquote("Meiotic genes in" ~
                 .(sub("_\\w+$", "", libName)) ~ "quantiles" ~
                 "(" * .(featRegion) * ") in" ~
                 .(gsub("_", " ", genomeName)) ~ .(region) ~
                 "(" * .(prettyNum(samplesNum,
                                   big.mark = ",",
                                   trim = T)) ~ " samples)"))
ggsave(paste0(plotDir,
              "bargraph_meiotic_gene_representation_among_", quantileLast,
              "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.pdf"),
       plot = bp,
       height = 8, width = 12)

