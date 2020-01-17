#!/applications/R/R-3.5.0/bin/Rscript

library(rtracklayer)
library(parallel)
library(data.table)
library(regioneR)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrStart <- rep(1, times = length(chrs))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"

NLRs <- readGFF(paste0(inDir, "NLRs_Krasileva/NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_representative_mRNA.gff3"))
NLRParser <- read.table(paste0(inDir, "NLRs_Krasileva/NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_representative_mRNA_NLRParser_out.tsv"),
                        header = T, stringsAsFactors = F)
NLRParser_complete <- NLRParser[NLRParser$complete == "complete",]
NLRParser_complete$SequenceName <- sub(pattern = "\\:\\:.+", replacement = "",
                                       x = NLRParser_complete$SequenceName)
NLRs_complete <- NLRs[NLRs$group %in% NLRParser_complete$SequenceName,]
NLRs_complete <- NLRs_complete[order(NLRs_complete$seqid,
                                    NLRs_complete$start,
                                    NLRs_complete$end),]
print(dim(NLRs_complete))
write.table(NLRs_complete[,1:9],
            file = paste0(inDir,
                          "NLRs_Krasileva/NB_ARC_genes_complete_NLRs_IWGSC_v1_Ksenia_Krasileva_representative_mRNA.gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
NLRs_complete_bed <- data.frame(chr = as.character(NLRs_complete[,1]),
                       start = as.integer(NLRs_complete[,4]-1),
                       end = as.integer(NLRs_complete[,5]),
                       name = as.character(NLRs_complete[,9]),
                       score = as.numeric(NLRs_complete[,6]),
                       strand = as.character(NLRs_complete[,7]),
                       stringsAsFactors = F)
write.table(NLRs_complete_bed,
            file = paste0(inDir,
                          "NLRs_Krasileva/NB_ARC_genes_complete_NLRs_IWGSC_v1_Ksenia_Krasileva_representative_mRNA.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
