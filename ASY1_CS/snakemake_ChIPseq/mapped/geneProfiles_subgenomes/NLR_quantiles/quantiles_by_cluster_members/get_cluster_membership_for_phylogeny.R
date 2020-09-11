#!/applications/R/R-3.5.0/bin/Rscript

# Get NLR physical cluster membership for each gene ID in the phylogeny
# (interactive Figure 3 of Steuernagel et al. 2020 Plant Physiol.)
# for determining whether physical clusters of NLRs are monophyletic

outDir <- "NLR_phylogeny_value_mapping/"
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

# Load features
features <- read.table("Asia/features_4quantiles_by_cluster_members_of_NLR_genes_in_Agenome_Bgenome_Dgenome_genomewide_Asia.txt",
                       header = T, sep = "\t", stringsAsFactors = F)
featureIDs <- sub(pattern = "\\.\\d+", replacement = "",
                  features$featureID)
features$cluster <- paste0(features$seqnames, "_c", features$cluster)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]

# Load representative gene model IDs as used for the RAxML phylogeny
# and replace our representative gene model IDs with these where they differ (89 genes)
NLR_phylo_IDs <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/NLRs_Steuernagel_Wulff_2020_Plant_Physiol/Figure_3_NLR_phylogeny_nexus_IDs.txt",
                            header = F, stringsAsFactors = F)[,1]
NLR_phylo_IDs <- NLR_phylo_IDs[!grepl(pattern = "TraesCSU", x = NLR_phylo_IDs)]
NLR_phylo_IDs <- NLR_phylo_IDs[!grepl(pattern = "LC.", x = NLR_phylo_IDs)]
NLR_phylo_IDs_alt_gene_models <- NLR_phylo_IDs[!(NLR_phylo_IDs %in% features$featureID)]
NLR_phylo_IDs_alt_gene_models <- NLR_phylo_IDs_alt_gene_models[which(substring(NLR_phylo_IDs_alt_gene_models, first = 1, last = 18) %in%
                                                                     substring(features$featureID, first = 1, last = 18))]

NLR_phylo_IDs_alt_gene_models_indices <- sapply(seq_along(NLR_phylo_IDs_alt_gene_models), function(x) {
     which(substring(NLR_phylo_IDs_alt_gene_models, first = 1, last = 18) %in%
           substring(as.character(features$featureID[which(substring(features$featureID, first = 1, last = 18) %in%
                                                           substring(NLR_phylo_IDs_alt_gene_models, first = 1, last = 18))])[x],
                     first = 1, last = 18))
})
NLR_phylo_IDs_alt_gene_models <- NLR_phylo_IDs_alt_gene_models[NLR_phylo_IDs_alt_gene_models_indices]
stopifnot(identical(substring(NLR_phylo_IDs_alt_gene_models, first = 1, last = 18),
                    substring(features$featureID[which(substring(features$featureID, first = 1, last = 18) %in%
                                                       substring(NLR_phylo_IDs_alt_gene_models, first = 1, last = 18))], first = 1, last = 18)))

features$featureID[which(substring(features$featureID, first = 1, last = 18) %in%
                         substring(NLR_phylo_IDs_alt_gene_models, first = 1, last = 18))] <- NLR_phylo_IDs_alt_gene_models

features <- features[features$featureID %in% NLR_phylo_IDs,]
features <- features[features$cluster_members > 1,]
# Create tables of physical cluster membership for phylogeny mapping
#cluster_membership_phylo_DF <- data.frame(Tree_node = paste0(sub(pattern = "\\.\\d+", replacement = "", x = features$featureID),
#                                                            "|", features$featureID),
#                                          Cluster = features$cluster,
#                                          stringsAsFactors = F)
cluster_membership_phylo_DF <- data.frame(Tree_node = features$featureID,
                                          Cluster = features$cluster,
                                          stringsAsFactors = F)
write.table(cluster_membership_phylo_DF,
            file = "Steuernagel_2020_Plant_Physiol_Figure_3_NLR_phylogeny_physical_clustering.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)
