#!/applications/R/R-3.5.0/bin/Rscript

# Write cM/Mb quantiles featureIDs for GO enrichment analysis

featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
                               split = ","))
quantiles <- 4

featuresDF <- read.table("quantiles_by_cMMb/WesternEurope/features_4quantiles_by_cMMb_of_genes_in_Agenome_Bgenome_Dgenome_genomewide_WesternEurope.txt",
                         header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)

# Get feature indices for each quantile
featureIndicesList <- lapply(seq_along(1:quantiles), function(k) {
  which(featuresDF$quantile == paste0("Quantile ", k))
})

featureIDsQuantileList <- lapply(seq_along(1:quantiles), function(k) {
  sub(pattern = "\\.\\d+", replacement = "",
      x = as.vector(featuresDF[featureIndicesList[[k]],]$featureID))
})
sapply(seq_along(featureIDsQuantileList), function(k) {
  write.table(featureIDsQuantileList[[k]],
              file = paste0("quantiles_by_cMMb/",
                            "featureIDs_quantile", k, "_of_", quantiles,
                            "_by_cMMb_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), ".txt"),
              quote = F, row.names = F, col.names = F)
})
