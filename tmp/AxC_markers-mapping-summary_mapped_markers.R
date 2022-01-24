#!/applications/R/R-4.0.0/bin/Rscript

# Compute signals of epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq,
# within physical marker intervals to enable correlations with wild type and fancm mutant
# crossover rates (cM/Mb)

# Usage:
# 

args <- commandArgs(trailingOnly = T)


options(stringsAsFactors = F)
library(data.table)

tab <- fread("AxC_markers-mapping-summary_mapped_markers.tsv")
colnames(tab) <- c("physical_marker", "chr", "pos", "wt_marker", "wt_cM", "fancm_marker", "fancm_cM")
tab$pos <- as.integer(gsub(",", "", tab$pos))
tab$chr <- sub("^", "chr", tab$chr)
  
