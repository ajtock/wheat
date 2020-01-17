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

mRNA_rep <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"))
print(dim(mRNA_rep))
#[1] 107891      9
mRNA_rep$groupTrunc <- sub(pattern = "\\.\\d+", replacement = "",
                           x = mRNA_rep$group)

NLR_IDs <- gsub(pattern = "^>", replacement = "",
                x = as.character(read.table(paste0(inDir, "NLRs_Krasileva/NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_v170120.ids"),
                                            header = F)$V1))
NLR_IDs <- gsub(pattern = "T----CS", replacement = "TraesCS",
                x = NLR_IDs)
NLRs <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_IDs,]
NLRs <- NLRs[order(NLRs$seqid,
                   NLRs$start,
                   NLRs$end),]
print(dim(NLRs))
write.table(NLRs[,1:9],
            file = paste0(inDir,
                          "NLRs_Krasileva/NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_representative_mRNA.gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
NLRs_bed <- data.frame(chr = as.character(NLRs[,1]),
                       start = as.integer(NLRs[,4]-1),
                       end = as.integer(NLRs[,5]),
                       name = as.character(NLRs[,9]),
                       score = as.numeric(NLRs[,6]),
                       strand = as.character(NLRs[,7]),
                       stringsAsFactors = F)
write.table(NLRs_bed,
            file = paste0(inDir,
                          "NLRs_Krasileva/NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_representative_mRNA.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
