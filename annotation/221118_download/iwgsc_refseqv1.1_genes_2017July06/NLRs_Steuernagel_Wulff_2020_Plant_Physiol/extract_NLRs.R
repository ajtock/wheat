#!/applications/R/R-3.5.0/bin/Rscript

# Extract NLR-encoding genes with matching IDs in IWGSC RefSeq v1.1 annotation
 
library(rtracklayer)

# Load NLR table
NLR <- read.table("Table_S5_NLRs_overlapping_refseq_v1.1_genes.txt",
                   header = T, sep = "\t", stringsAsFactors = F)
# Include NLR loci that overlap only 1 RefSeq v1.1 gene
NLR <- NLR[NLR$Number.of.overlapping.genes == 1,]  
NLR_IDs <- NLR$Overlapping.Genes
NLR_IDs <- unlist(strsplit(NLR_IDs, split = ","))
NLR_IDs <- sort(unique(NLR_IDs))
NLR_IDs_motifs <- sapply(NLR_IDs, function(x) {
  NLR[NLR$Overlapping.Genes == x,]$Motifs
})
NLR_IDs_motifs_multi_indices <- as.vector( which( 
  sapply(NLR_IDs_motifs, function(x) {
    length(x) > 1
  })
))
tmp <- sapply(NLR_IDs_motifs_multi_indices, function(x) {
  motifsSplit <- strsplit(NLR_IDs_motifs[[x]], split = ",")
  motifsSplitLens <- sapply(motifsSplit, function(y) length(y))
  motifsSplitMaxLen <- which(motifsSplitLens == max(motifsSplitLens))
  motifsSplitMaxLen[1]
})
sapply(seq_along(NLR_IDs_motifs), function(x) {
  print(length(NLR_IDs_motifs[[x]]))
})

# Get complete NLRs (with P-loop, NB-ARC and LRR domains)
NLR_complete <- NLR[grepl("complete", NLR$Completeness),]
NLR_complete_IDs <- NLR_complete$Overlapping.Genes
NLR_complete_IDs <- unlist(strsplit(NLR_complete_IDs, split = ","))
NLR_complete_IDs <- sort(unique(NLR_complete_IDs))
# Get partial NLRs
NLR_partial <- NLR[grepl("partial", NLR$Completeness),]
NLR_partial_IDs <- NLR_partial$Overlapping.Genes
NLR_partial_IDs <- unlist(strsplit(NLR_partial_IDs, split = ","))
NLR_partial_IDs <- sort(unique(NLR_partial_IDs))

# Load representative genes
inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
mRNA_rep <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"))
print(dim(mRNA_rep))
#[1] 107891      9
mRNA_rep$groupTrunc <- sub(pattern = "\\.\\d+", replacement = "",
                           x = mRNA_rep$group)

# Create tables of each set of NLR gene coordinates  in GFF3 and BED format
# all
NLR <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_IDs,]
NLR <- NLR[order(NLR$seqid,
                 NLR$start,
                 NLR$end),]
print(dim(NLR))
write.table(NLR[,1:9],
            file = "NLR_genes_representative_mRNA.gff3",
            quote = F, sep = "\t", row.names = F, col.names = F)
NLR_bed <- data.frame(chr = as.character(NLR[,1]),
                      start = as.integer(NLR[,4]-1),
                      end = as.integer(NLR[,5]),
                      name = as.character(NLR[,9]),
                      score = as.numeric(NLR[,6]),
                      strand = as.character(NLR[,7]),
                      stringsAsFactors = F)
write.table(NLR_bed,
            file = "NLR_genes_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)

# complete
NLR_complete <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_complete_IDs,]
NLR_complete <- NLR_complete[order(NLR_complete$seqid,
                                   NLR_complete$start,
                                   NLR_complete$end),]
print(dim(NLR_complete))
write.table(NLR_complete[,1:9],
            file = "NLR_genes_complete_representative_mRNA.gff3",
            quote = F, sep = "\t", row.names = F, col.names = F)
NLR_complete_bed <- data.frame(chr = as.character(NLR_complete[,1]),
                               start = as.integer(NLR_complete[,4]-1),
                               end = as.integer(NLR_complete[,5]),
                               name = as.character(NLR_complete[,9]),
                               score = as.numeric(NLR_complete[,6]),
                               strand = as.character(NLR_complete[,7]),
                               stringsAsFactors = F)
write.table(NLR_complete_bed,
            file = "NLR_genes_complete_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)

# partial
NLR_partial <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_partial_IDs,]
NLR_partial <- NLR_partial[order(NLR_partial$seqid,
                                 NLR_partial$start,
                                 NLR_partial$end),]
print(dim(NLR_partial))
write.table(NLR_partial[,1:9],
            file = "NLR_genes_partial_representative_mRNA.gff3",
            quote = F, sep = "\t", row.names = F, col.names = F)
NLR_partial_bed <- data.frame(chr = as.character(NLR_partial[,1]),
                              start = as.integer(NLR_partial[,4]-1),
                              end = as.integer(NLR_partial[,5]),
                              name = as.character(NLR_partial[,9]),
                              score = as.numeric(NLR_partial[,6]),
                              strand = as.character(NLR_partial[,7]),
                              stringsAsFactors = F)
write.table(NLR_partial_bed,
            file = "NLR_genes_partial_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)
