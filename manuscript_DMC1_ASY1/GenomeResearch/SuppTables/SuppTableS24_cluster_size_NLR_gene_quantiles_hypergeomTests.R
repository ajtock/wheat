#!/applications/R/R-3.5.0/bin/Rscript

# Compile NLR cluster size gene quantile hypergeometric test results into one supplemental table

library(parallel)

geneCat <- c(
             "NLR_PAMP_upreg_gene"
            )
geneCatName <- c(
                 "PAMP-triggered NLR-encoding"
                ) 
quantiles <- 1:4

summaryTableDFList <- mclapply(seq_along(geneCat), function(w) {
  summaryTableDFList_geneCat <- lapply(seq_along(quantiles), function(x) {
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/NLR_quantiles/",
                "quantiles_by_cluster_members/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_cluster_members_of_NLR_genes_in_Agenome_Bgenome_Dgenome_genomewide_hypergeomTestRes.RData"))
    Wg <- hgTestResults
    rm(hgTestResults); gc()
    Wg <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "All",
                     NLR_cluster_size_gene_quantile = paste0("Quantile ", x),
                     genes_in_quantile = Wg@quantile_genes,
                     expected = Wg@expected,
                     observed = Wg@observed,
                     alpha = Wg@alpha0.05,
                     log2obsexp = Wg@log2obsexp,
                     log2alphaexp = Wg@log2alpha,
                     pval = Wg@pval, 
                     stringsAsFactors = F)
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/NLR_quantiles/",
                "quantiles_by_cluster_members/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_cluster_members_of_NLR_genes_in_Agenome_genomewide_hypergeomTestRes.RData"))
    Ag <- hgTestResults
    rm(hgTestResults); gc()
    Ag <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "A",
                     NLR_cluster_size_gene_quantile = paste0("Quantile ", x),
                     genes_in_quantile = Ag@quantile_genes,
                     expected = Ag@expected,
                     observed = Ag@observed,
                     alpha = Ag@alpha0.05,
                     log2obsexp = Ag@log2obsexp,
                     log2alphaexp = Ag@log2alpha,
                     pval = Ag@pval, 
                     stringsAsFactors = F)
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/NLR_quantiles/",
                "quantiles_by_cluster_members/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_cluster_members_of_NLR_genes_in_Bgenome_genomewide_hypergeomTestRes.RData"))
    Bg <- hgTestResults
    rm(hgTestResults); gc()
    Bg <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "B",
                     NLR_cluster_size_gene_quantile = paste0("Quantile ", x),
                     genes_in_quantile = Bg@quantile_genes,
                     expected = Bg@expected,
                     observed = Bg@observed,
                     alpha = Bg@alpha0.05,
                     log2obsexp = Bg@log2obsexp,
                     log2alphaexp = Bg@log2alpha,
                     pval = Bg@pval, 
                     stringsAsFactors = F)
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/NLR_quantiles/",
                "quantiles_by_cluster_members/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_cluster_members_of_NLR_genes_in_Dgenome_genomewide_hypergeomTestRes.RData"))
    Dg <- hgTestResults
    rm(hgTestResults); gc()
    Dg <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "D",
                     NLR_cluster_size_gene_quantile = paste0("Quantile ", x),
                     genes_in_quantile = Dg@quantile_genes,
                     expected = Dg@expected,
                     observed = Dg@observed,
                     alpha = Dg@alpha0.05,
                     log2obsexp = Dg@log2obsexp,
                     log2alphaexp = Dg@log2alpha,
                     pval = Dg@pval, 
                     stringsAsFactors = F)
    rbind(Wg, Ag, Bg, Dg)
  })
  do.call(rbind, summaryTableDFList_geneCat)
}, mc.cores = length(geneCat), mc.preschedule = F)

summaryTableDF <- do.call(rbind, summaryTableDFList)

print(colnames(summaryTableDF))
colnames(summaryTableDF) <- c("NLR-encoding gene category", "Subgenome", "NLR cluster size gene quantile",
                              "Total NLR-encoding genes in quantile", "Expected category NLR-encoding genes in quantile", "Observed category NLR-encoding genes in quantile",
                              "Alpha (5%)", "Log2(observed/expected)", "Log2(alpha/expected)", "P")
print(colnames(summaryTableDF))
write.table(summaryTableDF,
            file = "Supplemental_TableS24_NLR_cluster_size_gene_quantiles_hypergeometricTests.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(summaryTableDF,
          file = "Supplemental_TableS24_NLR_cluster_size_gene_quantiles_hypergeometricTests.csv",
          row.names = F, quote = F)
