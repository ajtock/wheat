#!/applications/R/R-4.0.0/bin/Rscript

# Determine monophyly of physically clustered NLRs in the phylogeny
# (downloaded from iTOL interactive Figure 3 of Steuernagel et al. 2020 Plant Physiol.)

# Perform hypergeometric tests to determine whether each
# NLR-encoding gene quantile is over-represented or under-represented for
# NLRs that are members of conservatively-estimated monophyleyic physical clusters 
# (e.g., is the proportion of NLR genes within a given NLR gene quantile that
# are members of monophyleyic physical clusters significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_monophyletic) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_monophyletic) [m] + ( length(genome_genes) - length(genome_monophyletic)) [n]

# Usage
# ./analyse_monophyly_of_physically_clustered_NLRs_proportion_in_NLR_quantiles_hypergeometricTest.R 'cMMb' 'genes' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome' 100000

options(stringsAsFactors = FALSE)

#libName <- "cMMb"
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

if(libName %in% c("cMMb", "cluster_members")) {
  outDir <- paste0("../quantiles_by_", libName, "/hypergeometricTests/")
} else {
  outDir <- paste0("../quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", featRegion, "/hypergeometricTests/")
}
plotDir <- paste0(outDir, "plots/")

library(phylobase)
library(methods)
library(plotrix)
library(ggplot2)
library(ggbeeswarm)
#library(ggthemes)
library(grid)
library(gridExtra)
#library(extrafont)

# Load NLR phylogeny as phylo4 object
phylo <- readNewick(file = "Figure_3_NLR_phylogeny_newick.txt",
                    check.node.labels = "keep")
# Assign internal node numbers as labels (useful for various checks)
nodeLabels(phylo) <- paste("N", nodeId(phylo, "internal"), sep = "")

# Load cluster membership
clusters <- read.table("Steuernagel_2020_Plant_Physiol_Figure_3_NLR_phylogeny_physical_clustering.tsv",
                       header = T, sep = "\t")
# phlyobase needs all nodes to be present in clusters data
Tree_node <- c(clusters$Tree_node,
               as.character( phylo@label[which(!(phylo@label %in%
                                                 clusters$Tree_node))] ))
Cluster <- c(clusters$Cluster,
             rep(NA, length(as.character( phylo@label[which(!(phylo@label %in%
                                                                clusters$Tree_node))] ))))
clusters <- data.frame(Cluster = Cluster,
                       row.names = Tree_node)
clusterName <- unique(clusters$Cluster)[!is.na(unique(clusters$Cluster))]
## Order clustered NLR IDs so consistent with order in phylo@label
#label_indices <- sapply(seq_along(as.character(phylo@label)), function(x) {
#  which(clusters$Tree_node %in%
#        as.character(phylo@label)[x])
#})
#clusters <- clusters[label_indices,]

# Combine phylogeny with cluster membership
phylo_cl <- phylo4d(x = phylo,
                    all.data = clusters,
                    missing.data = "warn")
phylo_cl_df <- print(phylo_cl)

cluster_nodes <- lapply(seq_along(clusterName), function(x) {
  phylo_cl_df[which( !is.na(phylo_cl_df$Cluster) &
                            phylo_cl_df$Cluster == clusterName[x] ),]
})

cluster_nodes_MRCA_list <- list()
for(x in seq_along(cluster_nodes)) {
  if(dim(cluster_nodes[[x]])[1] > 1) {
    cluster_nodes_MRCA <- data.frame(label = cluster_nodes[[x]]$label, 
                                     node = cluster_nodes[[x]]$node,
                                     ancestor = cluster_nodes[[x]]$ancestor,
                                     MRCA = rep(MRCA(phylo, cluster_nodes[[x]]$node)[[1]],
                                                dim(cluster_nodes[[x]])[1]),
                                     MRCAlabel = rep(names(MRCA(phylo, cluster_nodes[[x]]$node)),
                                                     dim(cluster_nodes[[x]])[1]),
                                     cluster = cluster_nodes[[x]]$Cluster)
    cluster_nodes_MRCA_list <- c(cluster_nodes_MRCA_list, list(cluster_nodes_MRCA))
  }
}

cluster_nodes_MRCA_list_MRCA_among_ancestors <- lapply(seq_along(cluster_nodes_MRCA_list), function(x) {
  cluster_nodes_MRCA_list[[x]]$MRCA[1] %in% cluster_nodes_MRCA_list[[x]]$ancestor
})

cluster_nodes_MRCA_list[unlist(cluster_nodes_MRCA_list_MRCA_among_ancestors)]

print("Number of physical clusters of NLR genes whose most recent common ancestor is the parent node of at least one of those NLR genes within a given physical cluster =")
sum(unlist(cluster_nodes_MRCA_list_MRCA_among_ancestors))
#[1] 51
print("Number of physical clusters in NLR phylogeny with more than one NLR gene =")
length(cluster_nodes_MRCA_list_MRCA_among_ancestors)
#[1] 240
print("Proportion =")
sum(unlist(cluster_nodes_MRCA_list_MRCA_among_ancestors))/length(cluster_nodes_MRCA_list_MRCA_among_ancestors)
#[1] 0.2125

monophyletic_clusters <- do.call(rbind, cluster_nodes_MRCA_list[unlist(cluster_nodes_MRCA_list_MRCA_among_ancestors)])
monophyletic_clusters$label <- sub(pattern = "\\.\\d+", replacement = "",
                                   x = monophyletic_clusters$label)

# Perform hypergeometric test to analyse representation of
# monophyletic physically clustered NLRs among each NLR quantile
# Load feature quantiles
if(libName %in% c("cMMb", "cluster_members")) {
  featuresDF <- read.table(paste0(sub("hypergeometricTests/", "", outDir),
                                  "/WesternEurope/features_", quantileLast, "quantiles_by_",
                                  libName, "_of_NLR_genes_in_",
                                  genomeName, "_", region, "_WesternEurope.txt"),
                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
} else {
  featuresDF <- read.table(paste0(sub("hypergeometricTests/", "", outDir),
                                  "/WesternEurope/features_", quantileLast, "quantiles_by_",
                                  sub("_\\w+$", "", libName), "_in_", featRegion, "_of_NLR_genes_in_",
                                  genomeName, "_", region, "_WesternEurope.txt"),
                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)

}
featuresDF$featureID <- sub(pattern = "\\.\\d+", replacement = "",
                            x = featuresDF$featureID)
genome_genes <- featuresDF$featureID
quantile_genes_list <- lapply(quantileFirst:quantileLast, function(x) {
  featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", x) &
             featuresDF$featureID %in% genome_genes,]$featureID
})

featuresDF <- data.frame(featuresDF,
                         monophyletic = as.character(""))
for(x in 1:(dim(featuresDF)[1])) {
  if(featuresDF[x,]$featureID %in% monophyletic_clusters$label) {
    featuresDF[x,]$monophyletic <- "yes"
  } else {
    featuresDF[x,]$monophyletic <- "undetermined"
  }
}

# Get IDs of NLR genes that are members of monophyletic physical clusters
genome_monophyletic <- featuresDF[featuresDF$monophyletic == "yes",]$featureID
# Get the intersection of genome_monophyletic and genome_genes
# (this excludes genome_monophyletic genes not assigned to a chromosome)
genome_monophyletic <- intersect(genome_monophyletic, genome_genes)

quantile_genes_list <- lapply(quantileFirst:quantileLast, function(x) {
  featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", x) &
             featuresDF$featureID %in% genome_genes,]$featureID
})

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

# P-value is the probability of drawing >= length(quantile_monophyletic) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_monophyletic) [m] + ( length(genome_genes) - length(genome_monophyletic)) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(quantile_monophyletic) to length(quantile_genes)

#lapply(seq_along(quantile_genes_list), function(z) {
for(z in seq_along(quantile_genes_list)) {
  print(paste0("Quantile ", z))
  quantile_genes <- quantile_genes_list[[z]]
  # Get intersection of gene IDs in quantile z and gene IDs of NLRs
  quantile_monophyletic <- intersect(quantile_genes, genome_monophyletic)

  # Calculate the P-values for over-representation and under-representation
  # of NLRs among quantile z genes
  set.seed(2847502)
  # Over-representation:
  Pval_overrep <- sum(dhyper(x = length(quantile_monophyletic):length(quantile_genes),
                             m = length(genome_monophyletic),
                             n = length(genome_genes) - length(genome_monophyletic),
                             k = length(quantile_genes)))
  print(Pval_overrep)

  # Or by 1 minus the sum of the probabilities of drawing 0:(length(quantile_monophyletic)-1)
  print(1 - sum(dhyper(x = 0:(length(quantile_monophyletic)-1),
                       m = length(genome_monophyletic),
                       n = length(genome_genes) - length(genome_monophyletic),
                       k = length(quantile_genes))))

  # Under-representation
  Pval_underrep <- phyper(q = length(quantile_monophyletic),
                          m = length(genome_monophyletic),
                          n = length(genome_genes) - length(genome_monophyletic),
                          k = length(quantile_genes))
  print(Pval_underrep)

  # Sample without replacement
  hgDist <- rhyper(nn = samplesNum,
                   m = length(genome_monophyletic),
                   n = length(genome_genes) - length(genome_monophyletic),
                   k = length(quantile_genes))

  # Calculate P-values and significance levels
  if(length(quantile_monophyletic) > mean(hgDist)) {
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
                       observed = length(quantile_monophyletic),
                       expected = mean(hgDist),
                       log2obsexp = log2( (length(quantile_monophyletic)+1) / (mean(hgDist)+1) ),
                       log2alpha  = log2( (alpha0.05+1) / (mean(hgDist)+1) ),
                       quantile_genes = length(quantile_genes),
                       proportion_of_quantile = length(quantile_monophyletic) / length(quantile_genes),
                       random_proportions_of_quantile = hgDist / length(quantile_genes),
                       hypergeomDist = hgDist)
  if(libName %in% c("cMMb", "cluster_members")) {
  save(hgTestResults,
       file = paste0(outDir,
                     "NLR_clustered_monophyletic_gene_representation_among_quantile", z, "_of_", quantileLast,
                     "_by_", libName, "_of_NLR_genes_in_",
                     genomeName, "_", region, "_hypergeomTestRes.RData"))
  } else {
  save(hgTestResults,
       file = paste0(outDir,
                     "NLR_clustered_monophyletic_gene_representation_among_quantile", z, "_of_", quantileLast,
                     "_by_log2_", libName, "_control_in_", featRegion, "_of_NLR_genes_in_",
                     genomeName, "_", region, "_hypergeomTestRes.RData"))
  }

  # Generate histogram
  if(libName %in% c("cMMb", "cluster_members")) {
  pdf(paste0(plotDir,
             "NLR_clustered_monophyletic_gene_representation_among_quantile", z, "_of_", quantileLast,
             "_by_", libName, "_of_NLR_genes_in_",
             genomeName, "_", region, "_hypergeomTestRes_hist.pdf"),
             height = 4.5, width = 5)
  } else {
  pdf(paste0(plotDir,
             "NLR_clustered_monophyletic_gene_representation_among_quantile", z, "_of_", quantileLast,
             "_by_log2_", libName, "_control_in_", featRegion, "_of_NLR_genes_in_",
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
  titleText <- list(bquote("Monophyletic clustered NLRs in" ~
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
  if(libName %in% c("cMMb", "cluster_members")) {
  load(paste0(outDir,
              "NLR_clustered_monophyletic_gene_representation_among_quantile", z, "_of_", quantileLast,
              "_by_", libName, "_of_NLR_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.RData"))
  } else {
  load(paste0(outDir,
              "NLR_clustered_monophyletic_gene_representation_among_quantile", z, "_of_", quantileLast,
              "_by_log2_", libName, "_control_in_", featRegion, "_of_NLR_genes_in_",
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
  ggtitle(bquote("Monophyletic clustered NLRs in" ~
                 .(sub("_\\w+$", "", libName)) ~ "quantiles" ~
                 "(" * .(featRegion) * ") in" ~
                 .(gsub("_", " ", genomeName)) ~ .(region) ~
                 "(" * .(prettyNum(samplesNum,
                                   big.mark = ",",
                                   trim = T)) ~ " samples)"))
if(libName %in% c("cMMb", "cluster_members")) {
ggsave(paste0(plotDir,
              "bargraph_NLR_clustered_monophyletic_gene_representation_among_", quantileLast,
              "quantiles_by_", libName, "_of_NLR_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.pdf"),
       plot = bp,
       height = 8, width = 12)
} else {
ggsave(paste0(plotDir,
              "bargraph_NLR_clustered_monophyletic_gene_representation_among_", quantileLast,
              "quantiles_by_log2_", libName, "_control_in_", featRegion, "_of_NLR_genes_in_",
              genomeName, "_", region, "_hypergeomTestRes.pdf"),
       plot = bp,
       height = 8, width = 12)
}
