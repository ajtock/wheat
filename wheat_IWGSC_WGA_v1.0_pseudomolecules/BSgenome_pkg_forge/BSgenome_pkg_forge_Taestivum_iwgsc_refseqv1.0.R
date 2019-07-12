#!/applications/R/R-3.5.1/bin/Rscript

# Forge a BSgenome data package for the Triticum aestivum iwgsc_refseqv1.0 assembly
# ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/triticum_aestivum/dna/
# Consult https://bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf

# Usage:
# ./BSgenome_pkg_forge.R iwgsc_refseqv1.0

args <- commandArgs(trailingOnly = T)
# packageNamePt4 refers to the character string after the third "." in, e.g.,
# BSgenome.Taestivum.Cambridge.iwgsc_refseqv1.0 (i.e., iwgsc_refseqv1.0)
packageNamePt4 <- args[1]

library(BSgenome)

seed_files <- system.file("extdata", "HendersonLab", package = "BSgenome")
tail(list.files(seed_files, pattern = "-seed$"))

Ta_iwgsc_refseqv1.0_seed <- list.files(seed_files,
                                       pattern = paste0("\\.", packageNamePt4, "-seed$"),
                                       full.names = TRUE)
cat(readLines(Ta_iwgsc_refseqv1.0_seed), sep = "\n")

forgeBSgenomeDataPkg(Ta_iwgsc_refseqv1.0_seed)
