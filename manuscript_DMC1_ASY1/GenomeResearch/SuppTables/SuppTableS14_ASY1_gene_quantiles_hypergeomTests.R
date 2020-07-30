#!/applications/R/R-3.5.0/bin/Rscript

# Compile ASY1 gene quantile hypergeometric test results into one supplemental table

library(parallel)

geneCat <- c(
             "NLR_gene",
             "LAR_overlapping_gene",
             "MIR_overlapping_gene",
             "abiotic_stress_ASR_gene",
             "meiotic_gene",
             paste0("BIO0", 1:9, "_LAR_overlapping_gene"),
             paste0("BIO", 10:19, "_LAR_overlapping_gene")
            )
geneCatName <- c(
                 "NLR-encoding",
                 "LAR-overlapping",
                 "MIR-overlapping",
                 "Biotic and abiotic stress responsive",
                 "Meiotic",
                 "BIO1 (annual mean temp.) LAR-overlapping",
                 "BIO2 (mean diurnal range [mean of monthly (max. temp. - min. temp.)]) LAR-overlapping",
                 "BIO3 (isothermality [(BIO2/BIO7)x100]) LAR-overlapping",
                 "BIO4 (temp. seasonality (standard deviation x100) LAR-overlapping",
                 "BIO5 (max. temp. of warmest month) LAR-overlapping",
                 "BIO6 (min. temp. of coldest month) LAR-overlapping",
                 "BIO7 (temp. annual range [BIO5-BIO6]) LAR-overlapping",
                 "BIO8 (mean temp. of wettest quarter) LAR-overlapping",
                 "BIO9 (mean temp. of driest quarter) LAR-overlapping",
                 "BIO10 (mean temp. of warmest quarter) LAR-overlapping",
                 "BIO11 (mean temp. of coldest quarter) LAR-overlapping",
                 "BIO12 (annual precipitation) LAR-overlapping",
                 "BIO13 (precipitation of wettest month) LAR-overlapping",
                 "BIO14 (precipitation of driest month) LAR-overlapping",
                 "BIO15 (precipitation seasonality [coefficient of variation]) LAR-overlapping",
                 "BIO16 (precipitation of wettest quarter) LAR-overlapping",
                 "BIO17 (precipitation of driest quarter) LAR-overlapping",
                 "BIO18 (precipitation of warmest quarter) LAR-overlapping",
                 "BIO19 (precipitation of coldest quarter) LAR-overlapping"
                )
quantiles <- 1:4

summaryTableDFList <- mclapply(seq_along(geneCat), function(w) {
  summaryTableDFList_geneCat <- lapply(seq_along(quantiles), function(x) {
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                "quantiles_by_ASY1_in_genes/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_genes_of_genes_in_Agenome_Bgenome_Dgenome_genomewide_hypergeomTestRes.RData"))
    Wg <- hgTestResults
    rm(hgTestResults); gc()
    Wg <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "All",
                     ASY1_gene_quantile = paste0("Quantile ", x),
                     genes_in_quantile = Wg@quantile_genes,
                     expected = Wg@expected,
                     observed = Wg@observed,
                     alpha = Wg@alpha0.05,
                     log2obsexp = Wg@log2obsexp,
                     log2alphaexp = Wg@log2alpha,
                     pval = Wg@pval, 
                     stringsAsFactors = F)
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                "quantiles_by_ASY1_in_genes/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_genes_of_genes_in_Agenome_genomewide_hypergeomTestRes.RData"))
    Ag <- hgTestResults
    rm(hgTestResults); gc()
    Ag <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "A",
                     ASY1_gene_quantile = paste0("Quantile ", x),
                     genes_in_quantile = Ag@quantile_genes,
                     expected = Ag@expected,
                     observed = Ag@observed,
                     alpha = Ag@alpha0.05,
                     log2obsexp = Ag@log2obsexp,
                     log2alphaexp = Ag@log2alpha,
                     pval = Ag@pval, 
                     stringsAsFactors = F)
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                "quantiles_by_ASY1_in_genes/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_genes_of_genes_in_Bgenome_genomewide_hypergeomTestRes.RData"))
    Bg <- hgTestResults
    rm(hgTestResults); gc()
    Bg <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "B",
                     ASY1_gene_quantile = paste0("Quantile ", x),
                     genes_in_quantile = Bg@quantile_genes,
                     expected = Bg@expected,
                     observed = Bg@observed,
                     alpha = Bg@alpha0.05,
                     log2obsexp = Bg@log2obsexp,
                     log2alphaexp = Bg@log2alpha,
                     pval = Bg@pval, 
                     stringsAsFactors = F)
    load(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                "quantiles_by_ASY1_in_genes/hypergeometricTests/", geneCat[w], "_representation_among_quantile", x,
                "_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_genes_of_genes_in_Dgenome_genomewide_hypergeomTestRes.RData"))
    Dg <- hgTestResults
    rm(hgTestResults); gc()
    Dg <- data.frame(
                     gene_category = geneCatName[w],
                     subgenome = "D",
                     ASY1_gene_quantile = paste0("Quantile ", x),
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
colnames(summaryTableDF) <- c("Gene category", "Subgenome", "ASY1 gene quantile",
                              "Total genes in quantile", "Expected category genes in quantile", "Observed category genes in quantile",
                              "Alpha (5%)", "Log2(observed/expected)", "Log2(alpha/expected)", "P")
print(colnames(summaryTableDF))
write.table(summaryTableDF,
            file = "Supplemental_TableS14_ASY1_gene_quantiles_hypergeometricTests.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(summaryTableDF,
          file = "Supplemental_TableS14_ASY1_gene_quantiles_hypergeometricTests.csv",
          row.names = F, quote = F)
