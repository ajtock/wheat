#!/applications/R/R-3.5.0/bin/Rscript

# Compile DMC1 gene quantile hypergeometric test results into one supplemental table

library(parallel)

inDir <- "/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/quantiles_by_DMC1_in_genes/"
inDir1 <- paste0(inDir, "hypergeometricTests/homoeolog_exp_bias/HC_CS_no_stress/")
inDir2 <- paste0(inDir, "hypergeometricTests/triad_movement/HC_CS_no_stress/")

geneCat1 <- c("Balanced",
              ".dominant", "non_dominant", ".suppressed", "non_suppressed",
              "A.dominant", "B.dominant", "D.dominant",
              "A.suppressed", "B.suppressed", "D.suppressed")
geneCatName1 <- c("Balanced",
                  "Dominant", "Non-dominant", "Suppressed", "Non-suppressed",
                  "A dominant", "B dominant", "D dominant",
                  "A suppressed", "B suppressed", "D suppressed")
geneCat2 <- c("stable", "middle", "dynamic")
geneCatName2 <- c("Stable 10%", "Middle 80%", "Dynamic 10%")

geneCat <- c(geneCat1, geneCat2)
geneCatName <- c(geneCatName1, geneCatName2)

quantiles <- 1:4

summaryTableDFList <- mclapply(seq_along(quantiles), function(x) {
  summaryTableDFList_geneCat1 <- lapply(seq_along(geneCat1), function(w) {
    load(paste0(inDir1,
                "HC_CS_no_stress_", geneCat1[w], "_representation_among_quantile", x,
                "_of_4_by_log2_DMC1_Rep1_ChIP_control_in_genes_of_genes_in_Agenome_Bgenome_Dgenome_genomewide_hypergeomTestRes.RData"))
    WgCat1 <- hgTestResults
    rm(hgTestResults); gc()
    WgCat1 <- data.frame(
                         gene_category = geneCatName1[w],
                         DMC1_gene_quantile = paste0("Quantile ", x),
                         genes_in_quantile = WgCat1@quantile_genes,
                         expected = WgCat1@expected,
                         observed = WgCat1@observed,
                         alpha = WgCat1@alpha0.05,
                         log2obsexp = WgCat1@log2obsexp,
                         log2alphaexp = WgCat1@log2alpha,
                         pval = WgCat1@pval, 
                         stringsAsFactors = F)
    WgCat1
  })
  summaryTableDFList_geneCat2 <- lapply(seq_along(geneCat2), function(w) {
    load(paste0(inDir2,
                "HC_CS_no_stress_", geneCat2[w], "_representation_among_quantile", x,
                "_of_4_by_log2_DMC1_Rep1_ChIP_control_in_genes_of_genes_in_Agenome_Bgenome_Dgenome_genomewide_hypergeomTestRes_minConditions6.RData"))
    WgCat2 <- hgTestResults
    rm(hgTestResults); gc()
    WgCat2 <- data.frame(
                         gene_category = geneCatName2[w],
                         DMC1_gene_quantile = paste0("Quantile ", x),
                         genes_in_quantile = WgCat2@quantile_genes,
                         expected = WgCat2@expected,
                         observed = WgCat2@observed,
                         alpha = WgCat2@alpha0.05,
                         log2obsexp = WgCat2@log2obsexp,
                         log2alphaexp = WgCat2@log2alpha,
                         pval = WgCat2@pval, 
                         stringsAsFactors = F)
    WgCat2
  })
  rbind( do.call(rbind, summaryTableDFList_geneCat1),
         do.call(rbind, summaryTableDFList_geneCat2) )
}, mc.cores = length(quantiles), mc.preschedule = F)

summaryTableDF <- do.call(rbind, summaryTableDFList)

print(colnames(summaryTableDF))
colnames(summaryTableDF) <- c("Gene category", "DMC1 gene quantile",
                              "Total genes in quantile", "Expected category genes in quantile", "Observed category genes in quantile",
                              "Alpha (5%)", "Log2(observed/expected)", "Log2(alpha/expected)", "P")
print(colnames(summaryTableDF))
write.table(summaryTableDF,
            file = "Supplemental_TableS17_DMC1_gene_quantiles_hypergeometricTests_homoeologExpressionVariation.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(summaryTableDF,
          file = "Supplemental_TableS17_DMC1_gene_quantiles_hypergeometricTests_homoeologExpressionVariation.csv",
          row.names = F, quote = F)
