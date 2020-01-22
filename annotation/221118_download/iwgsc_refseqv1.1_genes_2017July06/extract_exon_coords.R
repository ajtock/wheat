#!/applications/R/R-3.5.0/bin/Rscript


# Create tables of exons within representative genes
# Write as gff and bed files

library(
library(rtracklayer)
library(parallel)
library(data.table)
library(regioneR)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = chrStart,
                                     end = chrLens),
                    strand = "*")
genomeGR <- genomeGR[grep(genomeName,
                          seqnames(genomeGR))@values]


chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrStart <- rep(1, times = length(chrs))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Load representative genes
inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
features <- readGFF(paste0(inDir,
                           "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"))
# Load genes including exons
# and get exons and introns for each representative gene model
genes <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706.gff3"))
genes <- genes[grep(genomeName, genes$seqid),]
exons <- genes[genes$type == "exon",]
exons <- exons[as.character(exons$Parent) %in% features$group,]
mclapply(seq_along(features$group), function(x) {
  exonParent <- exons[as.character(exons$Parent) == as.character(features$group)[x],]
}, mc.cores = detectCores())

mRNA <- genes[genes$type == "mRNA",]
print(dim(mRNA))
#[1] 43697    24


genes <- readGFF(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3"))
## Confirm that "first" exon start coordinate and last "exon" end coordinate
## equal mRNA start and end coordinates
#mclapply(seq_along(features$featureID), function(x) {
#  exonParent <- exons[as.character(exons$Parent) == as.character(features$featureID[x]),]
#  if( exonParent[1,]$start != features[features$featureID == features$featureID[x],]$start ) {
#    stop(paste0(as.character(features$featureID[x]),
#                ": first exon start coordinate does not equal mRNA start coordinate"))
#  }
#  if( exonParent[dim(exonParent)[1],]$end != features[features$featureID == features$featureID[x],]$end ) {
#    stop(paste0(as.character(features$featureID[x]),
#                ": last exon start coordinate does not equal mRNA end coordinate"))
#  }
#}, mc.cores = detectCores())
introns <- mclapply(seq_along(features$featureID), function(x) {
  exonParent <- exons[as.character(exons$Parent) == as.character(features$featureID[x]),]
  intronParent <- data.frame()
  for(i in 1:(dim(exonParent)[1]-1)) {
    exonParent[i,]$end


mRNA <- genes[genes$type == "mRNA",]
print(dim(mRNA))
#[1] 133744     24

# Obtain frequency of occurrence of each gene parent ID
n_occur <- data.frame(table(unlist(mRNA$Parent)))

# Obtain mRNA records for which the gene parent ID occurs only once
mRNA_unique <-  as.data.frame(
  mRNA[ unlist(mRNA$Parent)
    %in% n_occur$Var1[n_occur$Freq == 1],
  ]
)
# Obtain mRNA records for which the gene parent ID occurs more than once
mRNA_multi <- as.data.frame(
  mRNA[ unlist(mRNA$Parent)
    %in% n_occur$Var1[n_occur$Freq > 1],
  ]
)

# For each gene parent ID in mRNA_multi, obtain the mRNA record with the
# longest transcript
# If multiple mRNA records have the longest transcript,
# keep the first reported one only
mRNA_multi_list <- mclapply(seq_along(mRNA_multi[,1]), function(h) {
  mRNA_multi_ID_all <- mRNA_multi[ unlist(mRNA_multi$Parent)
                         == unlist(mRNA_multi[h,]$Parent),
                       ]
  mRNA_multi_ID_all[ mRNA_multi_ID_all$end-mRNA_multi_ID_all$start
    == max(mRNA_multi_ID_all$end-mRNA_multi_ID_all$start),
  ][1,]
}, mc.cores = detectCores())

# Collapse mRNA_multi_list into single data.frame and remove duplicates
mRNA_multi_dup <- rbindlist(mRNA_multi_list)
mRNA_multi_rep <- unique(as.data.frame(mRNA_multi_dup))

# Combine into one representative set of mRNA entries, order, and output as gff3
mRNA_rep <- rbind(mRNA_unique, mRNA_multi_rep)
mRNA_rep <- mRNA_rep[ order(mRNA_rep$seqid,
                            mRNA_rep$start,
                            mRNA_rep$end), ]
write.table(mRNA_rep[,1:9],
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_bed <- data.frame(chr = as.character(mRNA_rep[,1]),
                           start = as.integer(mRNA_rep[,4]-1),
                           end = as.integer(mRNA_rep[,5]),
                           name = as.integer(1:length(mRNA_rep[,1])),
                           score = as.numeric(mRNA_rep[,6]),
                           strand = as.character(mRNA_rep[,7]))
write.table(mRNA_rep_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_mRNA.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

mRNA_repGR <- GRanges(seqnames = mRNA_rep$seqid,
                      ranges = IRanges(start = mRNA_rep$start,
                                       end = mRNA_rep$end),
                      strand = mRNA_rep$strand)

# Generate coordinates for the same number of random loci
ranLocGR <- randomizeRegions(mRNA_repGR,
                             genome = genome,
                             per.chromosome = TRUE,
                             allow.overlaps = TRUE)
ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)))
write.table(ranLoc_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_mRNA_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
