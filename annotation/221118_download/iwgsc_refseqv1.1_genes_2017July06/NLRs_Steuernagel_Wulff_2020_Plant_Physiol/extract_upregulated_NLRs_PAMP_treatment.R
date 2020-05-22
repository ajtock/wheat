#!/applications/R/R-3.5.0/bin/Rscript

# Extract 266 NLR-encoding genes that were upregulated in response to PAMP treatment
# (i.e., those with a log2-fold-increase >= 1)

library(rtracklayer)

# Load gene expression table
geneexp <- read.csv("File_S1_gene_expression_PAMP_treatment.csv",
                    header = T, stringsAsFactors = F)
NLRexp <- geneexp[geneexp$nlr == "NLR",]
NLRexp <- NLRexp[,1:11]
NLR_ctn_upreg <- NLRexp[NLRexp$Wulff_C030 >= 1,]$gene
NLR_flg_upreg <- NLRexp[NLRexp$Wulff_F030 >= 1,]$gene
NLR_upreg_IDs <- sort(unique(c(NLR_ctn_upreg, NLR_flg_upreg)))

# Load representative genes
inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
mRNA_rep <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"))
print(dim(mRNA_rep))
#[1] 107891      9
mRNA_rep$groupTrunc <- sub(pattern = "\\.\\d+", replacement = "",
                           x = mRNA_rep$group)

NLR_upreg <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_upreg_IDs,]
NLR_upreg <- NLR_upreg[order(NLR_upreg$seqid,
                             NLR_upreg$start,
                             NLR_upreg$end),]
print(dim(NLR_upreg))
write.table(NLR_upreg[,1:9],
            file = "NLR_genes_upregulated_PAMP_treatment_representative_mRNA.gff3",
            quote = F, sep = "\t", row.names = F, col.names = F)
NLR_upreg_bed <- data.frame(chr = as.character(NLR_upreg[,1]),
                            start = as.integer(NLR_upreg[,4]-1),
                            end = as.integer(NLR_upreg[,5]),
                            name = as.character(NLR_upreg[,9]),
                            score = as.numeric(NLR_upreg[,6]),
                            strand = as.character(NLR_upreg[,7]),
                            stringsAsFactors = F)
write.table(NLR_upreg_bed,
            file = "NLR_genes_upregulated_PAMP_treatment_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)
