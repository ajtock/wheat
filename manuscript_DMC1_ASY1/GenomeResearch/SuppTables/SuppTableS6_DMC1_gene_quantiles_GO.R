#!/applications/R/R-3.5.0/bin/Rscript

# Compile DMC1 gene quantile GO enrichment analysis results into one supplemental table

library(GO.db)

quantiles <- 1:4

summaryTableDFList <- lapply(seq_along(quantiles), function(x) {
  Wg <- data.frame(Subgenome = "All",
                   DMC1_gene_quantile = paste0("Quantile ", x),
                   read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                                     "quantiles_by_log2_DMC1_Rep1_ChIP_control_in_genes/GO/featureIDs_quantile", x, "_of_4_by_log2_DMC1_Rep1_ChIP_control_in_genes_of_genes_in_",
                                     "Agenome_Bgenome_Dgenome_genomewide_GO_BP_enrichment.tsv"),
                                      header = T, sep = "\t", quote = "\"", check.names = F, stringsAsFactors = F),
                   stringsAsFactors = F)
  Ag <- data.frame(Subgenome = "A",
                   DMC1_gene_quantile = paste0("Quantile ", x),
                   read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                                     "quantiles_by_log2_DMC1_Rep1_ChIP_control_in_genes/GO/featureIDs_quantile", x, "_of_4_by_log2_DMC1_Rep1_ChIP_control_in_genes_of_genes_in_",
                                     "Agenome_genomewide_GO_BP_enrichment.tsv"),
                                      header = T, sep = "\t", quote = "\"", check.names = F, stringsAsFactors = F),
                   stringsAsFactors = F)
  Bg <- data.frame(Subgenome = "B",
                   DMC1_gene_quantile = paste0("Quantile ", x),
                   read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                                     "quantiles_by_log2_DMC1_Rep1_ChIP_control_in_genes/GO/featureIDs_quantile", x, "_of_4_by_log2_DMC1_Rep1_ChIP_control_in_genes_of_genes_in_",
                                     "Bgenome_genomewide_GO_BP_enrichment.tsv"),
                                      header = T, sep = "\t", quote = "\"", check.names = F, stringsAsFactors = F),
                   stringsAsFactors = F)
  Dg <- data.frame(Subgenome = "D",
                   DMC1_gene_quantile = paste0("Quantile ", x),
                   read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                                     "quantiles_by_log2_DMC1_Rep1_ChIP_control_in_genes/GO/featureIDs_quantile", x, "_of_4_by_log2_DMC1_Rep1_ChIP_control_in_genes_of_genes_in_",
                                     "Dgenome_genomewide_GO_BP_enrichment.tsv"),
                                      header = T, sep = "\t", quote = "\"", check.names = F, stringsAsFactors = F),
                   stringsAsFactors = F)
  rbind(Wg, Ag, Bg, Dg)
})
summaryTableDF <- do.call(rbind, summaryTableDFList)
# Remove additional test statistics from other methods ("classic", "elim", "weight")
summaryTableDF <- summaryTableDF[,-c(8:11)]
# Retain GO terms with topGO Fisher's exact test P-values <= 0.05
summaryTableDF$topGOFisher <- sub(pattern = "< ", replacement = "", summaryTableDF$topGOFisher)
summaryTableDF$topGOFisher <- as.numeric(summaryTableDF$topGOFisher)
summaryTableDF <- summaryTableDF[summaryTableDF$topGOFisher <= 0.05,]

# Replace truncated GO terms with full GO terms
summaryTableDF <- data.frame(summaryTableDF,
                             select(GO.db, keys = summaryTableDF$GO.ID,
                                    columns = c("GOID", "TERM"),#"DEFINITION"),
                                    keytype = "GOID"),
                             stringsAsFactors = F) 
# Ensure summaryTableDF$GO.ID matches the GO.db-obtained summaryTableDF$GOID,
# and thus truncated GO terms correspond to the correct full GO terms
stopifnot(all.equal(summaryTableDF$GO.ID, summaryTableDF$GOID))
print("Number of topGO-derived truncated GO terms to be replaced with corresponding GO.db-obtained full GO terms =")
print(length(summaryTableDF$Term[which(grepl("\\.\\.\\.", summaryTableDF$Term))]))
print("Number of non-matching summaryTableDF$Term (original) and summaryTableDF$TERM (GO.db-obtained) GO terms =")
print(all.equal(summaryTableDF$Term, summaryTableDF$TERM))
summaryTableDF$Term <- summaryTableDF$TERM
summaryTableDF <- summaryTableDF[,-c(9:10)]
summaryTableDF <- summaryTableDF[,c(1:5, 7, 6, 8)]
print(colnames(summaryTableDF))
colnames(summaryTableDF) <- c("Subgenome", "DMC1 gene quantile",
                              "GO ID", "GO term",
                              "Annotated genes in subgenome", "Expected genes in quantile", "Observed genes in quantile", "P")
print(colnames(summaryTableDF))
write.table(summaryTableDF,
            file = "Supplemental_TableS6_DMC1_gene_quantiles_GO.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(summaryTableDF,
          file = "Supplemental_TableS6_DMC1_gene_quantiles_GO.csv",
          row.names = F, quote = F)
