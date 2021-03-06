#!/applications/R/R-3.5.0/bin/Rscript

# Perform hypergeometric tests to determine whether each HudsonRM_all
# gene quantile is over-represented or under-represented for
# NLR-encoding genes (as defined by Ksenia Krasileva) and, separately,
# meiotic-module genes (as defined in Alabdullah et al. 2019 Front. Plant Sci.)
# (i.e., is the proportion of NLR-encoding genes or meiotic-module genes
# contained within each gene quantile significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_NLRs) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_NLRs) [m] + ( length(genome_genes) - length(genome_NLRs)) [n]

# Usage 
# ./proportion_NLRs_in_gene_quantiles_hypergeometricTest_HudsonRM.R 'HudsonRM_all' 'bodies' 1 2 'genomewide' 'Agenome_Bgenome_Dgenome' 100000

library(methods)
library(plotrix)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

#libName <- "HudsonRM_all"
#featRegion <- "bodies"
#quantileFirst <- 1
#quantileLast <- 2
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

outDir <- paste0("quantiles_by_", libName, "/hypergeometricTests/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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

for(p in 1:length(pop_name)) {
# Load feature quantiles
featuresDF <- read.table(paste0(sub("hypergeometricTests/", "", outDir),
                                "/", pop_name[p], "/features_", quantileLast, "quantiles_by_",
                                libName, "_of_genes_in_",
                                genomeName, "_", region, "_", pop_name[p], ".txt"),
                         header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
featuresDF$featureID <- sub(pattern = "\\.\\d+", replacement = "",
                            x = featuresDF$featureID)
genome_genes <- featuresDF$featureID
quantile_genes_list <- lapply(quantileFirst:quantileLast, function(x) {
  featuresDF[featuresDF$quantile == paste0("Quantile ", x),]$featureID
})
rm(featuresDF); gc()

# Load NLRs
#NLRs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
#                          "NLRs_Krasileva/NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_representative_mRNA.gff3"),
#                   header = F, stringsAsFactors = F)
NLRs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                          "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLR_genes_complete_representative_mRNA.gff3"),
                   header = F, stringsAsFactors = F)
chrs <- paste0(rep("chr", 21), rep(1:7, 3),
               c(rep("A", 7), rep("B", 7), rep("D", 7)))
genomeLetter <- unlist(strsplit(gsub("genome", "", genomeName), split = "_"))
# Subset NLRs to only those within a given subgenome
if(length(genomeLetter) == 1) {
  chrs <- chrs[grepl(genomeLetter, chrs)]
  NLRs <- NLRs[NLRs$V1 %in% chrs,]
}

# Replace gene model ID decimal suffix (e.g., ".1")
NLRs$V9 <- sub(pattern = "\\.\\d+", replacement = "",
               x = NLRs$V9)
genome_NLRs <- as.character(NLRs$V9)

genome_NLRs <- intersect(genome_NLRs, genome_genes)

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

# P-value is the probability of drawing >= length(quantile_NLRs) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_NLRs) [m] + ( length(genome_genes) - length(genome_NLRs)) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(quantile_NLRs) to length(quantile_genes)

#lapply(seq_along(quantile_genes_list), function(z) {
for(z in seq_along(quantile_genes_list)) {
  quantile_genes <- quantile_genes_list[[z]]
  # Get intersection of gene IDs in quantile z and gene IDs of NLRs
  quantile_NLRs <- intersect(quantile_genes, genome_NLRs)

  # Calculate the P-values for over-representation and under-representation
  # of NLRs among quantile z genes
  set.seed(2847502)
  # Over-representation:
  Pval_overrep <- sum(dhyper(x = length(quantile_NLRs):length(quantile_genes),
                             m = length(genome_NLRs),
                             n = length(genome_genes) - length(genome_NLRs),
                             k = length(quantile_genes)))
  print(Pval_overrep)

  # Or by 1 minus the sum of the probabilities of drawing 0:(length(quantile_NLRs)-1)
  print(1 - sum(dhyper(x = 0:(length(quantile_NLRs)-1),
                       m = length(genome_NLRs),
                       n = length(genome_genes) - length(genome_NLRs),
                       k = length(quantile_genes))))

  # Under-representation
  Pval_underrep <- phyper(q = length(quantile_NLRs),
                          m = length(genome_NLRs),
                          n = length(genome_genes) - length(genome_NLRs),
                          k = length(quantile_genes))
  print(Pval_underrep)

  # Sample without replacement
  hgDist <- rhyper(nn = samplesNum,
                   m = length(genome_NLRs),
                   n = length(genome_genes) - length(genome_NLRs),
                   k = length(quantile_genes))

  # Calculate P-values and significance levels
  if(length(quantile_NLRs) > mean(hgDist)) {
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
                       observed = length(quantile_NLRs),
                       expected = mean(hgDist),
                       log2obsexp = log2( length(quantile_NLRs) / mean(hgDist) ),
                       log2alpha  = log2( alpha0.05 / mean(hgDist) ),
                       quantile_genes = length(quantile_genes),
                       proportion_of_quantile = length(quantile_NLRs) / length(quantile_genes),
                       random_proportions_of_quantile = hgDist / length(quantile_genes),
                       hypergeomDist = hgDist)
  save(hgTestResults,
       file = paste0(outDir,
                     "NLR_gene_representation_among_quantile", z, "_of_", quantileLast,
                     "_by_", libName, "_of_genes_in_",
                     genomeName, "_", region, "_hypergeomTestRes_", pop_name[p], ".RData"))

  # Generate histogram
  pdf(paste0(plotDir,
             "NLR_gene_representation_among_quantile", z, "_of_", quantileLast,
             "_by_", libName, "_of_genes_in_",
             genomeName, "_", region, "_hypergeomTestRes_", pop_name[p], "_hist.pdf"),
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
        text = bquote("Genes (" * .(pop_name_plot[p]) * ")"),
        line = 1.85)
  titleText <- list(bquote("NLR-encoding genes in" ~
                           .(libName) ~ "Quantile" ~ .(as.character(z)) ~
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
  load(paste0(outDir,
              "NLR_gene_representation_among_quantile", z, "_of_", quantileLast,
              "_by_", libName, "_of_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes_", pop_name[p], ".RData"))
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
  ggtitle(bquote("NLR-encoding genes in" ~
                 .(libName) ~ "quantiles" ~
                 "(" * .(pop_name_plot[p]) * ") in" ~
                 .(gsub("_", " ", genomeName)) ~ .(region) ~
                 "(" * .(prettyNum(samplesNum,
                                   big.mark = ",",
                                   trim = T)) ~ " samples)"))
ggsave(paste0(plotDir,
              "bargraph_NLR_gene_representation_among_", quantileLast,
              "quantiles_by_", libName, "_of_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes_", pop_name[p], ".pdf"),
       plot = bp,
       height = 8, width = 12)
}
