#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./get_GO_IDs.R genomewide "IWGSC+Stress"

#region <- "genomewide"
#annoType <- "IWGSC+Stress"

args <- commandArgs(trailingOnly = T)
region <- args[1]
annoType <- args[2]

library(parallel)
library(dplyr)
library(tidyr)
library(tidyverse)

# Redundant given conditional statement below
e <- "Variable annoType is not one of 'GO', 'PO', 'TO', 'andrea_go', 'BUSCO', 'IWGSC+Stress', 'slim_GO', 'slim_andrea_go', 'slim_BUSCO', 'slim_IWGSC+Stress'"
tryCatch(stopifnot(annoType %in% c('GO', 'PO', 'TO', 'andrea_go', 'BUSCO', 'IWGSC+Stress', 'slim_GO', 'slim_andrea_go', 'slim_BUSCO', 'slim_IWGSC+Stress')),
         error = function(e) stop(e))

# Load table of functional annotations
anno <- readRDS("OntologiesForGenes.rds")

# Retain rows corresponding to annoType
anno <- anno[anno$ontology == annoType,]
anno <- anno[,1:2]

# Replace "1G" with "2G" in gene IDs for consistency with v1.1
anno$Gene <- sub(pattern = "1G", replacement = "2G",
                 x = anno$Gene)
# Spread so that each column is a gene ID, within which GO terms
# are listed (1 per row)
# See post by dmi3k at https://community.rstudio.com/t/spread-why-errors/2076/7
tmp1 <- anno %>%
  dplyr::group_by_at(vars(-ID)) %>% # group by everything other than the value column
  dplyr::mutate(row_id = 1:n()) %>% ungroup() %>% # build group index
  tidyr::spread(key = Gene, value = ID) %>% # spread
  dplyr::select(-row_id) # drop the index

tmp2 <- as.data.frame(t(tmp1))
tmp3 <- tidyr::unite(data = tmp2, col = "annoID",
                     sep = ",", remove = T)
anno <- data.frame(geneID = rownames(tmp3),
                   annoID = gsub(",NA", "", tmp3$annoID),
                   stringsAsFactors = F)

# Load representative genes in GFF3 format
genesGFF <- as.character(read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA.gff3")$V9)
genesGFF <- sub(pattern = "\\.\\d+", replacement = "",
                x = genesGFF)
# Retain annotations only for those gene IDs present in genesGFF,
# and add gene IDs present in genesGFF but absent from anno
# for fair enrichment analysis
anno <- anno[anno$geneID %in% genesGFF,]
absentIDs <- genesGFF[!(genesGFF %in% anno$geneID)]
absentIDsDF <- data.frame(geneID = absentIDs,
                          annoID = "")
anno <- rbind(anno,
              absentIDsDF)
anno <- anno[order(anno$geneID, decreasing = F),]

# Remove rows for gene IDs assigned to unanchored scaffolds
anno <- anno[!grepl("CSU", anno$geneID),]
write.table(anno,
            file = paste0("RamirezGonzalez_2018_iwgsc_refseqv1.0_OntologiesForGenes_FunctionalAnnotation_HCgenes_in_Agenome_Bgenome_Dgenome_",
                          region, "_GO_IDs_no_chrUn.tsv"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# Subset by subgenome and region
# A
AgenesGFF <- as.character(read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_Agenome_",
                                            region, ".gff3"))$V9)
AgenesGFF <- sub(pattern = "\\.\\d+", replacement = "",
                 x = AgenesGFF)
annoA <- anno[anno$geneID %in% AgenesGFF,]
write.table(annoA,
            file = paste0("RamirezGonzalez_2018_iwgsc_refseqv1.0_OntologiesForGenes_FunctionalAnnotation_HCgenes_in_Agenome_",
                          region, "_GO_IDs_no_chrUn.tsv"),
            row.names = F, col.names = F, quote = F, sep = "\t")
# B
BgenesGFF <- as.character(read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_Bgenome_",
                                            region, ".gff3"))$V9)
BgenesGFF <- sub(pattern = "\\.\\d+", replacement = "",
                 x = BgenesGFF)
annoB <- anno[anno$geneID %in% BgenesGFF,]
write.table(annoB,
            file = paste0("RamirezGonzalez_2018_iwgsc_refseqv1.0_OntologiesForGenes_FunctionalAnnotation_HCgenes_in_Bgenome_",
                          region, "_GO_IDs_no_chrUn.tsv"),
            row.names = F, col.names = F, quote = F, sep = "\t")
# D
DgenesGFF <- as.character(read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_Dgenome_",
                                            region, ".gff3"))$V9)
DgenesGFF <- sub(pattern = "\\.\\d+", replacement = "",
                 x = DgenesGFF)
annoD <- anno[anno$geneID %in% DgenesGFF,]
write.table(annoD,
            file = paste0("RamirezGonzalez_2018_iwgsc_refseqv1.0_OntologiesForGenes_FunctionalAnnotation_HCgenes_in_Dgenome_",
                          region, "_GO_IDs_no_chrUn.tsv"),
            row.names = F, col.names = F, quote = F, sep = "\t")
