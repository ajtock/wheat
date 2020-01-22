#!/applications/R/R-3.5.0/bin/Rscript

# Create tables of exons and introns within representative genes
# for each wheat subgenome and generate random loci of the same
# number and width distribution
# Write as GFF3 and BED files

# Usage:
# ./extract_exon_intron_coords_subgenomes.R 'genomewide' 'A'

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

# Load representative genes
inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
mRNA_rep <- readGFF(paste0(inDir,
                           "IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                           genomeName, "genome_", region, ".gff3"))
# Load genes including exons
# and get exons and introns for each representative gene model
genes <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706.gff3"))
genes <- genes[grep(genomeName, genes$seqid),]
exons <- genes[genes$type == "exon",]
exons_rep <- exons[as.character(exons$Parent) %in% mRNA_rep$group,]

# NOTE: when analysis region-specific mRNAs, must remove those in masked regions
# as done below for exons
# Remove exons in masked regions
exons_repGR <- GRanges(seqnames = exons_rep$seqid,
                       ranges = IRanges(start = exons_rep$start,
                                        end = exons_rep$end),
                       strand = exons_rep$strand,
                       source = exons_rep$source,
                       type = exons_rep$type,
                       score = exons_rep$score,
                       phase = exons_rep$phase,
                       ID = exons_rep$ID,
                       Parent = exons_rep$Parent)
exons_mask_overlaps <- findOverlaps(query = maskGR,
                                    subject = exons_repGR,
                                    ignore.strand = T,
                                    select = "all")
if(length(exons_mask_overlaps) > 0) {
  exons_repGR <- exons_repGR[-subjectHits(exons_mask_overlaps)]
}
exons_rep_gff <- data.frame(seqid = as.character(seqnames(exons_repGR)),
                            source = as.character(exons_repGR$source),
                            type = as.character(exons_repGR$type),
                            start = as.integer(start(exons_repGR)),
                            end = as.integer(end(exons_repGR)),
                            score = as.numeric(exons_repGR$score),
                            strand = as.character(strand(exons_repGR)),
                            phase = as.integer(exons_repGR$phase),
                            ID = as.character(exons_repGR$ID),
                            Parent = as.character(exons_repGR$Parent),
                            stringsAsFactors = F)
write.table(exons_rep_gff,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_exons_in_",
                          genomeName, "genome_", region, ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
                       
exons_rep_bed <- data.frame(chr = as.character(exons_rep_gff[,1]),
                            start = as.integer(exons_rep_gff[,4]-1),
                            end = as.integer(exons_rep_gff[,5]),
                            name = as.integer(1:length(exons_rep_gff[,1])),
                            score = as.numeric(exons_rep_gff[,6]),
                            strand = as.character(exons_rep_gff[,7]),
                            stringsAsFactors = F)
write.table(exons_rep_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_exons_in_",
                          genomeName, "genome_", region, ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Define function to select randomly positioned loci of the same
# width distribution as exons_repGR
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
# ranLocGR contains the same number of loci per chromosome as exons_repGR
chrs <- chrs[grep(genomeName, chrs)]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  exons_repChrGR <- exons_repGR[seqnames(exons_repGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(exons_repChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = (start(regionChrGR):end(regionChrGR)),
                                      n = length(exons_repChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(exons_repChrGR)),
                         strand = strand(exons_repChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)),
                         stringsAsFactors = F)
write.table(ranLoc_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_exons_in_",
                          genomeName, "genome_", region, "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Get introns
introns_rep_gff_list <- mclapply(seq_along(as.character(mRNA_rep$group)), function(x) {
#introns_rep_gff <- data.frame()
#for(x in seq_along(as.character(mRNA_rep$group))) {
#  print(x)
  exonParent <- exons_rep_gff[as.character(exons_rep_gff$Parent) == as.character(mRNA_rep$group)[x],]
  if(dim(exonParent)[1] > 1) {
    intronParent <- data.frame()
    for(i in 1:(dim(exonParent)[1]-1)) {
#      print(i)
      intron_i_start <- exonParent[i,]$end+1
      intron_i_end <- exonParent[i+1,]$start-1
      intron_i <- exonParent[i,]
      intron_i$type <- "intron"
      intron_i$ID <- gsub(pattern = "exon", replacement = "intron",
                          x = intron_i$ID)
      intron_i$start <- intron_i_start
      intron_i$end <- intron_i_end
      intronParent <- rbind(intronParent, intron_i)
    }
#    introns_rep_gff <- rbind(introns_rep_gff, intronParent)
    intronParent
  }
#}
}, mc.cores = detectCores())
introns_rep_gff <- do.call(rbind, introns_rep_gff_list)
write.table(introns_rep_gff,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_introns_in_",
                          genomeName, "genome_", region, ".gff3"),
            quote = F, sep = "\t", row.names = F, col.names = F)
                       
introns_rep_bed <- data.frame(chr = as.character(introns_rep_gff[,1]),
                              start = as.integer(introns_rep_gff[,4]-1),
                              end = as.integer(introns_rep_gff[,5]),
                              name = as.integer(1:length(introns_rep_gff[,1])),
                              score = as.numeric(introns_rep_gff[,6]),
                              strand = as.character(introns_rep_gff[,7]),
                              stringsAsFactors = F)
write.table(introns_rep_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_introns_in_",
                          genomeName, "genome_", region, ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

introns_repGR <- GRanges(seqnames = introns_rep_gff$seqid,
                         ranges = IRanges(start = introns_rep_gff$start,
                                          end = introns_rep_gff$end),
                         strand = introns_rep_gff$strand,
                         source = introns_rep_gff$source,
                         type = introns_rep_gff$type,
                         score = introns_rep_gff$score,
                         phase = introns_rep_gff$phase,
                         ID = introns_rep_gff$ID,
                         Parent = introns_rep_gff$Parent)

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as introns_repGR
chrs <- chrs[grep(genomeName, chrs)]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  introns_repChrGR <- introns_repGR[seqnames(introns_repGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(introns_repChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = (start(regionChrGR):end(regionChrGR)),
                                      n = length(introns_repChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(introns_repChrGR)),
                         strand = strand(introns_repChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)),
                         stringsAsFactors = F)
write.table(ranLoc_bed,
            file = paste0(inDir,
                          "IWGSC_v1.1_HC_20170706_representative_introns_in_",
                          genomeName, "genome_", region, "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

### Confirm that "first" exon start coordinate and last "exon" end coordinate
### equal mRNA start and end coordinates
##mclapply(seq_along(mRNA_rep$group), function(x) {
##  exonParent <- exons[as.character(exons$Parent) == as.character(mRNA_rep$group[x]),]
##  if( exonParent[1,]$start != mRNA_rep[mRNA_rep$group == mRNA_rep$group[x],]$start ) {
##    stop(paste0(as.character(mRNA_rep$group[x]),
##                ": first exon start coordinate does not equal mRNA start coordinate"))
##  }
##  if( exonParent[dim(exonParent)[1],]$end != mRNA_rep[mRNA_rep$group == mRNA_rep$group[x],]$end ) {
##    stop(paste0(as.character(mRNA_rep$group[x]),
##                ": last exon start coordinate does not equal mRNA end coordinate"))
##  }
##}, mc.cores = detectCores())
