#!/applications/R/R-3.5.0/bin/Rscript

# For all three wheat subgenomes, load and plot as bargraphs results from
# hypergeometric tests to determine whether each HudsonRM_all
# gene quantile is over-represented or under-represented for
# a given feature category; e.g.,
# NLR-encoding genes (as defined by Ksenia Krasileva) and, separately,
# meiotic-module genes (as defined in Alabdullah et al. 2019 Front. Plant Sci.)
# (i.e., is the propotion of NLR-encoding genes or meiotic-module genes
# contained within each gene quantile significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_NLRs) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_NLRs) [m] + ( length(genome_genes) - length(genome_NLRs)) [n]

# Usage 
# ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only_HudsonRM.R 'HudsonRM_all' 'bodies' 1 2 'genomewide' NLR 'NLR-encoding' 100000 'black,navy,dodgerblue4,deepskyblue'
# ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only_HudsonRM.R 'HudsonRM_all' 'bodies' 1 2 'genomewide' meiotic 'Meiotic' 100000 'black,purple4,purple,magenta'
# ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only_HudsonRM.R 'HudsonRM_all' 'bodies' 1 2 'genomewide' LAR_overlapping 'LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
# ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only_HudsonRM.R 'HudsonRM_all' 'bodies' 1 2 'genomewide' MIR_overlapping 'MIR-overlapping' 100000 'black,red4,red,tomato'

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
#featCat <- "NLR_overlapping"
#featCatPlot <- "NLR-overlapping"
#samplesNum <- 100000
#genomeColours <- unlist(strsplit('black,darkgreen,seagreen,springgreen', split = ",")) 

args <- commandArgs(trailingOnly = TRUE)
libName <- args[1]
featRegion <- args[2]
quantileFirst <- as.integer(args[3])
quantileLast <- as.integer(args[4])
region <- args[5]
featCat <- args[6]
featCatPlot <- args[7]
samplesNum <- as.numeric(args[8])
genomeColours <- unlist(strsplit(args[9], split = ","))

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

options(scipen = 100)

# Plot bar graph summarising permutation test results
genomeNames <- c("Agenome_Bgenome_Dgenome", "Agenome", "Bgenome", "Dgenome")
for(p in 1:length(pop_name)) {
hg_list <- lapply(seq_along(genomeNames), function(y) {
  hg_list_quantile <- list() 
  for(z in quantileFirst:quantileLast) {
    load(paste0(outDir,
                featCat, "_gene_representation_among_quantile", z, "_of_", quantileLast,
                "_by_", libName, "_of_genes_in_", genomeNames[y], "_",
                region, "_hypergeomTestRes_", pop_name[p], ".RData"))
    hg_list_quantile <- c(hg_list_quantile, hgTestResults)
  }
  return(hg_list_quantile)
})
bargraph_df <- data.frame(Subgenome = rep(c("All genomes", "A genome", "B genome", "D genome"), each = quantileLast),
                          Quantile = rep(paste0("Quantile ", quantileFirst:quantileLast), 4),
                          log2ObsExp = c(sapply(seq_along(genomeNames), function(y) {
                                           sapply(seq_along(hg_list[[y]]), function(x) {
                                             hg_list[[y]][[x]]@log2obsexp
                                           })
                                         })),
                          log2alpha0.05 = c(sapply(seq_along(genomeNames), function(y) {
                                              sapply(seq_along(hg_list[[y]]), function(x) {
                                                hg_list[[y]][[x]]@log2alpha
                                              })
                                            })))
bargraph_df$Quantile <- factor(bargraph_df$Quantile,
                               levels = paste0("Quantile ", quantileFirst:quantileLast))
bargraph_df$Subgenome <- factor(bargraph_df$Subgenome,
                                levels = c("All genomes", "A genome", "B genome", "D genome"))
bp <- ggplot(data = bargraph_df,
             mapping = aes(x = Quantile,
                           y = log2ObsExp,
                           fill = Subgenome)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = genomeColours,
                    labels = levels(bargraph_df$Subgenome)) +
  geom_point(mapping = aes(x = Quantile,
                           y = log2alpha0.05),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey80", size = 20) +
  labs(y = bquote("Log"[2]*"(observed/expected) genes in quantile")) +
  scale_y_continuous(limits = c(-0.6, 0.3)) +
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
        #legend.position = "none",
        #legend.position = c(0.05, 0.30),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 40.5, 5.5), "pt"),
        plot.title = element_text(size = 18, colour = "black", hjust = 0.5)) +
  ggtitle(bquote(.(featCatPlot) ~ "genes in" ~
                 .(libName) ~ "quantiles" ~
                 "(" * .(pop_name_plot[p]) * ")" ~
                 "(" * .(prettyNum(samplesNum,
                                   big.mark = ",",
                                   trim = T)) ~ "samples)"))
ggsave(paste0(plotDir,
              "bargraph_", featCat, "_gene_representation_among_", quantileLast,
              "quantiles_by_", libName, "_of_genes_in_each_subgenome_",
              region, "_hypergeomTestRes_", pop_name[p], ".pdf"),
       plot = bp,
       height = 8, width = 16)
}
