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
# Get motifs annotated to each NLR
NLR_IDs_motifs <- sapply(NLR_IDs, function(x) {
  NLR[NLR$Overlapping.Genes == x,]$Motifs
})
# Get indices of NLRs annotated with more than 1 set of motifs
# This occurs because some genes in the RefSeq v1.1 annotation
# are overlapped by more than one NLR in the Steuernagel et al. (2020) Plant Physiol. analysis,
# each containing a potentially different set of motifs
NLR_IDs_motifs_multi_indices <- as.vector( which( 
  sapply(NLR_IDs_motifs, function(x) {
    length(x) > 1
  })
))
# If the aim is to evaluate representation of the number of motifs
# among NLRs in different NLR quantiles (e.g., grouped by decreasing cM/Mb),
# then it is necessary to identify a representative set of motifs (i.e., the largest set)
# for each gene, to avoid double-counting the same motif occurrence represented in two sets
# Alternatively, if the aim is to evaluate representation of NLRs that contain a given motif
# among NLRs in different NLR quantiles, then it would be preferable to concatenate motif sets
NLR_IDs_repres_motifs <- NLR_IDs_motifs
NLR_IDs_concat_motifs <- NLR_IDs_motifs 
# For each of the above indices, get the index of the longest set of motifs 
NLR_IDs_motifs_multi_maxLen_indices <- sapply(NLR_IDs_motifs_multi_indices, function(x) {
  motifsSplit <- strsplit(NLR_IDs_motifs[[x]], split = ",")
  motifsSplitLens <- sapply(motifsSplit, function(y) length(y))
  motifsSplitMaxLen <- which(motifsSplitLens == max(motifsSplitLens))
  motifsSplitMaxLen[1]
})
# For each of the NLRs annotated with more than 1 set of motifs,
# assign the longest set of motifs
for(x in seq_along(NLR_IDs_motifs_multi_indices)) {
  NLR_IDs_repres_motifs[[NLR_IDs_motifs_multi_indices[x]]] <- NLR_IDs_repres_motifs[[NLR_IDs_motifs_multi_indices[x]]][NLR_IDs_motifs_multi_maxLen_indices[x]]
}
# For each of the NLRs annotated with more than 1 set of motifs,
# concatenate the sets
for(x in seq_along(NLR_IDs_motifs_multi_indices)) {
  NLR_IDs_concat_motifs[[NLR_IDs_motifs_multi_indices[x]]] <- paste0(NLR_IDs_concat_motifs[[NLR_IDs_motifs_multi_indices[x]]], collapse = ",")
}

# To enable easy pattern matching, flank motif numbers with "_"
for(x in seq_along(NLR_IDs_repres_motifs)) {
  subunderx <- gsub(pattern = ",", replacement = "_",
                    x = NLR_IDs_repres_motifs[x])
  subunderx <- gsub("^(\\d)", "_\\1", subunderx)
  subunderx <- gsub("(\\d)$", "\\1_", subunderx)
  NLR_IDs_repres_motifs[[x]] <- subunderx
}
for(x in seq_along(NLR_IDs_concat_motifs)) {
  subunderx <- gsub(pattern = ",", replacement = "_",
                    x = NLR_IDs_concat_motifs[x])
  subunderx <- gsub("^(\\d)", "_\\1", subunderx)
  subunderx <- gsub("(\\d)$", "\\1_", subunderx)
  NLR_IDs_concat_motifs[[x]] <- subunderx
}

 
# Get complete NLRs (with P-loop, NB-ARC and LRR domains)
NLR_complete <- NLR[grepl("complete", NLR$Completeness),]
NLR_complete_IDs <- NLR_complete$Overlapping.Genes
NLR_complete_IDs <- unlist(strsplit(NLR_complete_IDs, split = ","))
NLR_complete_IDs <- sort(unique(NLR_complete_IDs))
# Get motifs annotated to each NLR
NLR_complete_IDs_motifs <- sapply(NLR_complete_IDs, function(x) {
  NLR[NLR$Overlapping.Genes == x,]$Motifs
})
# Get indices of NLRs annotated with more than 1 set of motifs
# This occurs because some genes in the RefSeq v1.1 annotation
# are overlapped by more than one NLR in the Steuernagel et al. (2020) Plant Physiol. analysis,
# each containing a potentially different set of motifs
NLR_complete_IDs_motifs_multi_indices <- as.vector( which( 
  sapply(NLR_complete_IDs_motifs, function(x) {
    length(x) > 1
  })
))
# If the aim is to evaluate representation of the number of motifs
# among NLRs in different NLR quantiles (e.g., grouped by decreasing cM/Mb),
# then it is necessary to identify a representative set of motifs (i.e., the largest set)
# for each gene, to avoid double-counting the same motif occurrence represented in two sets
# Alternatively, if the aim is to evaluate representation of NLRs that contain a given motif
# among NLRs in different NLR quantiles, then it would be preferable to concatenate motif sets
NLR_complete_IDs_repres_motifs <- NLR_complete_IDs_motifs
NLR_complete_IDs_concat_motifs <- NLR_complete_IDs_motifs 
# For each of the above indices, get the index of the longest set of motifs 
NLR_complete_IDs_motifs_multi_maxLen_indices <- sapply(NLR_complete_IDs_motifs_multi_indices, function(x) {
  motifsSplit <- strsplit(NLR_complete_IDs_motifs[[x]], split = ",")
  motifsSplitLens <- sapply(motifsSplit, function(y) length(y))
  motifsSplitMaxLen <- which(motifsSplitLens == max(motifsSplitLens))
  motifsSplitMaxLen[1]
})
# For each of the NLRs annotated with more than 1 set of motifs,
# assign the longest set of motifs
for(x in seq_along(NLR_complete_IDs_motifs_multi_indices)) {
  NLR_complete_IDs_repres_motifs[[NLR_complete_IDs_motifs_multi_indices[x]]] <- NLR_complete_IDs_repres_motifs[[NLR_complete_IDs_motifs_multi_indices[x]]][NLR_complete_IDs_motifs_multi_maxLen_indices[x]]
}
# For each of the NLRs annotated with more than 1 set of motifs,
# concatenate the sets
for(x in seq_along(NLR_complete_IDs_motifs_multi_indices)) {
  NLR_complete_IDs_concat_motifs[[NLR_complete_IDs_motifs_multi_indices[x]]] <- paste0(NLR_complete_IDs_concat_motifs[[NLR_complete_IDs_motifs_multi_indices[x]]], collapse = ",")
}

# To enable easy pattern matching, flank motif numbers with "_"
for(x in seq_along(NLR_complete_IDs_repres_motifs)) {
  subunderx <- gsub(pattern = ",", replacement = "_",
                    x = NLR_complete_IDs_repres_motifs[x])
  subunderx <- gsub("^(\\d)", "_\\1", subunderx)
  subunderx <- gsub("(\\d)$", "\\1_", subunderx)
  NLR_complete_IDs_repres_motifs[[x]] <- subunderx
}
for(x in seq_along(NLR_complete_IDs_concat_motifs)) {
  subunderx <- gsub(pattern = ",", replacement = "_",
                    x = NLR_complete_IDs_concat_motifs[x])
  subunderx <- gsub("^(\\d)", "_\\1", subunderx)
  subunderx <- gsub("(\\d)$", "\\1_", subunderx)
  NLR_complete_IDs_concat_motifs[[x]] <- subunderx
}

# Get partial NLRs
NLR_partial <- NLR[grepl("partial", NLR$Completeness),]
NLR_partial_IDs <- NLR_partial$Overlapping.Genes
NLR_partial_IDs <- unlist(strsplit(NLR_partial_IDs, split = ","))
NLR_partial_IDs <- sort(unique(NLR_partial_IDs))
# Get motifs annotated to each NLR
NLR_partial_IDs_motifs <- sapply(NLR_partial_IDs, function(x) {
  NLR[NLR$Overlapping.Genes == x,]$Motifs
})
# Get indices of NLRs annotated with more than 1 set of motifs
# This occurs because some genes in the RefSeq v1.1 annotation
# are overlapped by more than one NLR in the Steuernagel et al. (2020) Plant Physiol. analysis,
# each containing a potentially different set of motifs
NLR_partial_IDs_motifs_multi_indices <- as.vector( which( 
  sapply(NLR_partial_IDs_motifs, function(x) {
    length(x) > 1
  })
))
# If the aim is to evaluate representation of the number of motifs
# among NLRs in different NLR quantiles (e.g., grouped by decreasing cM/Mb),
# then it is necessary to identify a representative set of motifs (i.e., the largest set)
# for each gene, to avoid double-counting the same motif occurrence represented in two sets
# Alternatively, if the aim is to evaluate representation of NLRs that contain a given motif
# among NLRs in different NLR quantiles, then it would be preferable to concatenate motif sets
NLR_partial_IDs_repres_motifs <- NLR_partial_IDs_motifs
NLR_partial_IDs_concat_motifs <- NLR_partial_IDs_motifs 
# For each of the above indices, get the index of the longest set of motifs 
NLR_partial_IDs_motifs_multi_maxLen_indices <- sapply(NLR_partial_IDs_motifs_multi_indices, function(x) {
  motifsSplit <- strsplit(NLR_partial_IDs_motifs[[x]], split = ",")
  motifsSplitLens <- sapply(motifsSplit, function(y) length(y))
  motifsSplitMaxLen <- which(motifsSplitLens == max(motifsSplitLens))
  motifsSplitMaxLen[1]
})
# For each of the NLRs annotated with more than 1 set of motifs,
# assign the longest set of motifs
for(x in seq_along(NLR_partial_IDs_motifs_multi_indices)) {
  NLR_partial_IDs_repres_motifs[[NLR_partial_IDs_motifs_multi_indices[x]]] <- NLR_partial_IDs_repres_motifs[[NLR_partial_IDs_motifs_multi_indices[x]]][NLR_partial_IDs_motifs_multi_maxLen_indices[x]]
}
# For each of the NLRs annotated with more than 1 set of motifs,
# concatenate the sets
for(x in seq_along(NLR_partial_IDs_motifs_multi_indices)) {
  NLR_partial_IDs_concat_motifs[[NLR_partial_IDs_motifs_multi_indices[x]]] <- paste0(NLR_partial_IDs_concat_motifs[[NLR_partial_IDs_motifs_multi_indices[x]]], collapse = ",")
}

# To enable easy pattern matching, flank motif numbers with "_"
for(x in seq_along(NLR_partial_IDs_repres_motifs)) {
  subunderx <- gsub(pattern = ",", replacement = "_",
                    x = NLR_partial_IDs_repres_motifs[x])
  subunderx <- gsub("^(\\d)", "_\\1", subunderx)
  subunderx <- gsub("(\\d)$", "\\1_", subunderx)
  NLR_partial_IDs_repres_motifs[[x]] <- subunderx
}
for(x in seq_along(NLR_partial_IDs_concat_motifs)) {
  subunderx <- gsub(pattern = ",", replacement = "_",
                    x = NLR_partial_IDs_concat_motifs[x])
  subunderx <- gsub("^(\\d)", "_\\1", subunderx)
  subunderx <- gsub("(\\d)$", "\\1_", subunderx)
  NLR_partial_IDs_concat_motifs[[x]] <- subunderx
}


# Load representative genes
inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
mRNA_rep <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"))
print(dim(mRNA_rep))
#[1] 107891      9
mRNA_rep$groupTrunc <- sub(pattern = "\\.\\d+", replacement = "",
                           x = mRNA_rep$group)

# Create tables of each set of NLR gene coordinates in GFF3 and BED format
# all
NLR <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_IDs,]
NLR <- NLR[order(NLR$groupTrunc),]
NLR_IDs_repres_motifs_subset <- NLR_IDs_repres_motifs[names(NLR_IDs_repres_motifs) %in%
                                                            NLR$groupTrunc]
NLR_IDs_repres_motifs_subset <- NLR_IDs_repres_motifs_subset[order(names(NLR_IDs_repres_motifs_subset))]
stopifnot(identical(NLR$groupTrunc, names(NLR_IDs_repres_motifs_subset)))
stopifnot(all.equal(NLR$groupTrunc, names(NLR_IDs_repres_motifs_subset)))
NLR$score <- as.character(NLR_IDs_repres_motifs_subset)
                  
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
                      score = as.character(NLR[,6]),
                      strand = as.character(NLR[,7]),
                      stringsAsFactors = F)
write.table(NLR_bed,
            file = "NLR_genes_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)

# complete
NLR_complete <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_complete_IDs,]
NLR_complete <- NLR_complete[order(NLR_complete$groupTrunc),]
NLR_complete_IDs_repres_motifs_subset <- NLR_complete_IDs_repres_motifs[names(NLR_complete_IDs_repres_motifs) %in%
                                                                        NLR_complete$groupTrunc]
NLR_complete_IDs_repres_motifs_subset <- NLR_complete_IDs_repres_motifs_subset[order(names(NLR_complete_IDs_repres_motifs_subset))]
stopifnot(identical(NLR_complete$groupTrunc, names(NLR_complete_IDs_repres_motifs_subset)))
stopifnot(all.equal(NLR_complete$groupTrunc, names(NLR_complete_IDs_repres_motifs_subset)))
NLR_complete$score <- as.character(NLR_complete_IDs_repres_motifs_subset)
                  
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
                               score = as.character(NLR_complete[,6]),
                               strand = as.character(NLR_complete[,7]),
                               stringsAsFactors = F)
write.table(NLR_complete_bed,
            file = "NLR_genes_complete_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)

# partial
NLR_partial <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_partial_IDs,]
NLR_partial <- NLR_partial[order(NLR_partial$groupTrunc),]
NLR_partial_IDs_repres_motifs_subset <- NLR_partial_IDs_repres_motifs[names(NLR_partial_IDs_repres_motifs) %in%
                                                                      NLR_partial$groupTrunc]
NLR_partial_IDs_repres_motifs_subset <- NLR_partial_IDs_repres_motifs_subset[order(names(NLR_partial_IDs_repres_motifs_subset))]
stopifnot(identical(NLR_partial$groupTrunc, names(NLR_partial_IDs_repres_motifs_subset)))
stopifnot(all.equal(NLR_partial$groupTrunc, names(NLR_partial_IDs_repres_motifs_subset)))
NLR_partial$score <- as.character(NLR_partial_IDs_repres_motifs_subset)
                  
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
                              score = as.character(NLR_partial[,6]),
                              strand = as.character(NLR_partial[,7]),
                              stringsAsFactors = F)
write.table(NLR_partial_bed,
            file = "NLR_genes_partial_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)
