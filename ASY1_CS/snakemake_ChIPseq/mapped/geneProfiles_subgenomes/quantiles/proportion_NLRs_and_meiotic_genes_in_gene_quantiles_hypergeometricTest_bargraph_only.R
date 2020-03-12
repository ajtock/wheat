#!/applications/R/R-3.5.0/bin/Rscript

# Perform hypergeometric tests to determine whether each ASY1 or DMC1
# gene quantile is over-represented or under-represented for
# NLR-encoding genes (as defined by Ksenia Krasileva) and, separately,
# meiotic-module genes (as defined in Alabdullah et al. 2019 Front. Plant Sci.)
# (i.e., is the propotion of NLR-encoding genes or meiotic-module genes
# contained within each gene quantile significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_NLRs) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_NLRs) [m] + ( length(genome_genes) - length(genome_NLRs)) [n]

### Perform hypergeometric test to determine whether a significant proportion of 
### elements belonging to each transposon family are differentially expressed in kyp suvh5 suvh6

# P-value is the probability of drawing >= length(DE_TEfamIDs) [x] TEs
# in a sample size of length(DE_TEIDs) [k] from a total TE set consisting of
# length(TAIR10_TEfamIDs) [m] + (length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs)) [n]

# Usage 
# ./proportion_NLRs_and_meiotic_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R 'ASY1_CS_Rep1_ChIP' 'bodies' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome' 100000

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

options(scipen = 100)

# Plot bar graph summarising permutation test results
pt_list_NLRs <- list()
for(z in quantileFirst:quantileLast) {
  load(
       paste0(outDir,
              "NLR_gene_representation_among_quantile", z, "_of_", quantileLast,
              "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.RData"))
  pt_list_NLRs <- c(pt_list_NLRs, hgTestResults)
}
pt_list_meio <- list()
for(z in quantileFirst:quantileLast) {
  load(
       paste0(outDir,
              "meiotic_gene_representation_among_quantile", z, "_of_", quantileLast,
              "_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.RData"))
  pt_list_meio <- c(pt_list_meio, hgTestResults)
}
bargraph_df <- data.frame(Feature = rep(c("NLR-encoding genes", "Meiotic genes"), each = quantileLast),
                          Quantile = rep(paste0("Quantile ", quantileFirst:quantileLast), 2),
                          log2ObsExp = c(sapply(seq_along(pt_list_NLRs), function(x) { pt_list_NLRs[[x]]@log2obsexp }),
                                         sapply(seq_along(pt_list_meio), function(x) { pt_list_meio[[x]]@log2obsexp })),
                          log2alpha0.05 = c(sapply(seq_along(pt_list_NLRs), function(x) { pt_list_NLRs[[x]]@log2alpha }),
                                            sapply(seq_along(pt_list_meio), function(x) { pt_list_meio[[x]]@log2alpha })))
bargraph_df$Quantile <- factor(bargraph_df$Quantile,
                               levels = paste0("Quantile ", quantileFirst:quantileLast))
bargraph_df$Feature <- factor(bargraph_df$Feature,
                              levels = c("NLR-encoding genes", "Meiotic genes"))
bp <- ggplot(data = bargraph_df,
             mapping = aes(x = Quantile,
                           y = log2ObsExp,
                           fill = Feature)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("dodgerblue",
                               "darkorange"),
                    labels = levels(bargraph_df$Feature)) +
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
  ggtitle(bquote("NLR-encoding and meiotic genes in" ~
                 .(sub("_\\w+$", "", libName)) ~ "quantiles" ~
                 "(" * .(featRegion) * ") in" ~
                 .(gsub("_", " ", genomeName)) ~ .(region) ~
                 "(" * .(prettyNum(samplesNum,
                                   big.mark = ",",
                                   trim = T)) ~ " samples)"))
ggsave(paste0(plotDir,
              "bargraph_NLR_and_meiotic_gene_representation_among_", quantileLast,
              "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.pdf"),
       plot = bp,
       height = 8, width = 12)

