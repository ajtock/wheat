#!/applications/R/R-3.5.0/bin/Rscript

# Create tables of representative genes for each wheat subgenome
# and generate random loci of the same number and width distribution
# Write as gff and bed files 

# Usage:
# ./extract_mRNA_coords_subgenomes.R 'genomewide' 'A'

args <- commandArgs(trailingOnly = T)
region <- args[1]
genomeName <- args[2]

library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(data.table)

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

# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "pericentromeres") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "pericentromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genome
wide")
}

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"

genes <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706.gff3"))
genes <- genes[grep(genomeName, genes$seqid),]
mRNA <- genes[genes$type == "mRNA",]
print(dim(mRNA))
#[1] 43697    24

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

# Combine into one representative set of mRNA entries, order,
# and output in GFF3 and BED formats
mRNA_rep <- rbind(mRNA_unique, mRNA_multi_rep)
mRNA_rep <- mRNA_rep[ order(mRNA_rep$seqid,
                            mRNA_rep$start,
                            mRNA_rep$end), ]
#write.table(mRNA_rep[,1:9],
#            file = paste0(inDir,
#                          "IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
#                          genomeName, "genome.gff3"),
#            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_bed <- data.frame(chr = as.character(mRNA_rep[,1]),
                           start = as.integer(mRNA_rep[,4]-1),
                           end = as.integer(mRNA_rep[,5]),
                           name = as.integer(1:length(mRNA_rep[,1])),
                           score = as.numeric(mRNA_rep[,6]),
                           strand = as.character(mRNA_rep[,7]))
write.table(mRNA_rep_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                          genomeName, "genome_", region, ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

mRNA_repGR <- GRanges(seqnames = mRNA_rep$seqid,
                      ranges = IRanges(start = mRNA_rep$start,
                                       end = mRNA_rep$end),
                      strand = mRNA_rep$strand)

# Define function to select randomly positioned loci of the same
# width distribution as mRNA_repGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as mRNA_repGR
chrs <- chrs[grep(genomeName, chrs)]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  mRNA_repChrGR <- mRNA_repGR[seqnames(mRNA_repGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(mRNA_repChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = (start(regionChrGR):end(regionChrGR)),
                                      n = length(mRNA_repChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(mRNA_repChrGR)),
                         strand = strand(mRNA_repChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)))
write.table(ranLoc_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_mRNA_randomLoci_in_",
                          genomeName, "genome_", region, ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
