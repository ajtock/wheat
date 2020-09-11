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
# ./analyse_monophyly_of_physically_clustered_NLRs.R 'cMMb' 'genes' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome' 100000

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

sum(unlist(cluster_nodes_MRCA_list_MRCA_among_ancestors))
#[1] 51
length(cluster_nodes_MRCA_list_MRCA_among_ancestors)
#[1] 240
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

# P-value is the probability of drawing >= length(quantile_monophyletic) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_monophyletic) [m] + ( length(genome_genes) - length(genome_monophyletic)) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(quantile_monophyletic) to length(quantile_genes)

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
  print("P-value for over-representation =")
  print(Pval_overrep)

  # Or by 1 minus the sum of the probabilities of drawing 0:(length(quantile_monophyletic)-1)
  print("P-value for over-representation =")
  print(1 - sum(dhyper(x = 0:(length(quantile_monophyletic)-1),
                       m = length(genome_monophyletic),
                       n = length(genome_genes) - length(genome_monophyletic),
                       k = length(quantile_genes))))

  # Under-representation
  print("P-value for under-representation =")
  Pval_underrep <- phyper(q = length(quantile_monophyletic),
                          m = length(genome_monophyletic),
                          n = length(genome_genes) - length(genome_monophyletic),
                          k = length(quantile_genes))
  print(Pval_underrep)
}



#
#
#
#
#
#


