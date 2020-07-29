#!/applications/R/R-3.5.0/bin/Rscript

# Compile ASY1 peak overlap permutation test results into one supplemental table

# Compile alignment stats for multiple data sets into one supplemental table

subgenomes <- c("A", "B", "D")

summaryTableDFList <- lapply(seq_along(subgenomes), function(x) {
  R1R3 <- data.frame(Subgenome_compartments = paste0(subgenomes[x], "-genome R1 & R3"),
                     rbind(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/",
                                             "PeakRanger1.18/ranger/p0.001_q0.01/euchromatin/regioneR/noMinWidth_mergedOverlaps/",
                                             "permTest_10000perms_ASY1_CS_Rep1_ChIP_peaks_vs_others_in_",
                                             subgenomes[x],"genome_euchromatin_DataFrame.txt"), header = T, stringsAsFactors = F),
                           read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/",
                                             "PeakRanger1.18/ranger/p0.001_q0.01/euchromatin/regioneR/noMinWidth_mergedOverlaps/",
                                             "permTest_10000perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
                                             subgenomes[x],"genome_euchromatin_DataFrame.txt"), header = T, stringsAsFactors = F)),
                     stringsAsFactors = F)
  R2C <- data.frame(Subgenome_compartments = paste0(subgenomes[x], "-genome R2a-R2b"),
                    rbind(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/",
                                            "PeakRanger1.18/ranger/p0.001_q0.01/heterochromatin/regioneR/noMinWidth_mergedOverlaps/",
                                            "permTest_10000perms_ASY1_CS_Rep1_ChIP_peaks_vs_others_in_",
                                            subgenomes[x],"genome_heterochromatin_DataFrame.txt"), header = T, stringsAsFactors = F),
                          read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/",
                                            "PeakRanger1.18/ranger/p0.001_q0.01/heterochromatin/regioneR/noMinWidth_mergedOverlaps/",
                                            "permTest_10000perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
                                            subgenomes[x],"genome_heterochromatin_DataFrame.txt"), header = T, stringsAsFactors = F)),
                    stringsAsFactors = F)
  rbind(R1R3, R2C)
})
summaryTableDF <- do.call(rbind, summaryTableDFList)

# Amend feature names for inclusion in table
featureNames <- summaryTableDF$feature
featureNames <- gsub(pattern = "promoters", replacement = "gene_promoters", x = featureNames)
featureNames <- gsub(pattern = "terminators", replacement = "gene_terminators", x = featureNames)
featureNames <- gsub(pattern = "TSSsPlus500bp", replacement = "gene_5'_ends", x = featureNames)
featureNames <- gsub(pattern = "TTSsMinus500bp", replacement = "gene_3'_ends", x = featureNames)
featureNames <- gsub(pattern = "NLRgene", replacement = "NLR-encoding_gene", x = featureNames)
featureNames <- gsub(pattern = "meiogene", replacement = "Meiotic_gene", x = featureNames)
featureNames <- gsub(pattern = "^NLRs$", replacement = "NLR-encoding_genes", x = featureNames)
featureNames <- gsub(pattern = "^meio$", replacement = "Meiotic_genes", x = featureNames)
featureNames <- gsub(pattern = "^gene", replacement = "Gene", x = featureNames)
featureNames <- gsub(pattern = "H2AZ", replacement = "H2A.Z", x = featureNames)
featureNames <- sub(pattern = "_\\w\\w\\w$", replacement = "", x = featureNames)
featureNames <- gsub(pattern = "_", replacement = " ", x = featureNames)
summaryTableDF$feature <- featureNames
colnames(summaryTableDF) <- c("Subgenome compartments", "Feature", "Feature number",
                              "Expected overlaps", "Observed overlaps",
                              "P", "Local Z-score")
write.table(summaryTableDF,
            file = "Supplemental_TableS3_ASY1_peak_overlap_summary.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(summaryTableDF,
          file = "Supplemental_TableS3_ASY1_peak_overlap_summary.csv",
          row.names = F, quote = F)
