#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 28.04.2020

# Separate peaks into those in the following wheat subgenome compartments (a.k.a. partitions):
# 1. R1 and R3 (distal)
# 2. R2a and R2b (interstitial)
# 3. C (proximal)
# Write in GFF and BED formats 

# Usage:
# ./get_peaks_in_compartments.R ASY1_CS_Rep1_ChIP 'distal' 'A'

#libName <- "ASY1_CS_Rep1_ChIP"
#region <- "distal"
#genomeName <- "A"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
region <- args[2]
genomeName <- args[3]

library(GenomicRanges)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                               "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,2])
centromereEnd <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                             "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = chrStart,
                                     end = chrLens),
                    strand = "*")
genomeGR <- genomeGR[grep(genomeName,
                          seqnames(genomeGR))@values]

# Define region to be analysed
if(region == "distal") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "interstitial") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(chrPartitions$R1_R2a+1,
                                                 chrPartitions$C_R2b),
                                       end = c(chrPartitions$R2a_C,
                                               chrPartitions$R2b_R3-1)),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "proximal") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
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
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else {
  stop("region is not distal, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "distal") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "interstitial") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 3),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2a_C+1,
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrPartitions$C_R2b-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "proximal") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
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
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else {
  stop("region is not distal, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Import peaks as GRanges object
load(paste0(libName, "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
peaksGR <- peaksGR[grep(genomeName,
                        seqnames(peaksGR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaksGR <- peaksGR[-subjectHits(mask_peaks_overlap)]
strand(peaksGR) <- "*"
print("***********peaks***********")
print(peaksGR)

# Define function to select randomly positioned loci of the same
# width distribution as peaksGR
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
# ranLocGR contains the same number of loci per chromosome as peaksGR
chrs <- chrs[grep(genomeName, chrs)]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  peaksChrGR <- peaksGR[seqnames(peaksGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(peaksChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {
                                                             start(regionChrGR[x]) : end(regionChrGR[x])
                                                           })),
                                      n = length(peaksChrGR))
  ranLocChrIR <- IRanges(start = ranLocChrStart,
                         width = width(peaksChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = ranLocChrIR,
                         strand = "*")
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}

# Convert into GFF3 and BED formats
# peaks
peaks <- data.frame(peaksGR)

peaksgff <- data.frame(chr = as.character(peaks$seqnames),
                       source = as.character(rep(".")),
                       feature = as.character(rep(paste0(libName, "_", region, "_peak"))),
                       start = as.integer(peaks$start),
                       end = as.integer(peaks$end),
                       score = as.character(rep(".")),
                       strand = as.character(rep(".")),
                       frame = as.character(rep(".")),
                       attribute = as.character(rep(".")))
write.table(peaksgff,
            file = paste0(libName,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, ".gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

peaksbed <- data.frame(chr = as.character(peaks$seqnames),
                       start = as.integer(peaks$start-1),
                       end = as.integer(peaks$end),
                       name = as.integer(1:length(peaks$seqnames)),
                       score = rep("NA", length(peaks$seqnames)),
                       strand = as.character(peaks$strand))
write.table(peaksbed,
            file = paste0(libName,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, ".bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")


# ranLoc
ranLoc <- data.frame(ranLocGR)

ranLocgff <- data.frame(chr = as.character(ranLoc$seqnames),
                        source = as.character(rep(".")),
                        feature = as.character(rep(paste0(libName, "_", region, "_ranLoc"))),
                        start = as.integer(ranLoc$start),
                        end = as.integer(ranLoc$end),
                        score = as.character(rep(".")),
                        strand = as.character(rep(".")),
                        frame = as.character(rep(".")),
                        attribute = as.character(rep(".")))
write.table(ranLocgff,
            file = paste0(libName,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, "_randomLoci.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

ranLocbed <- data.frame(chr = as.character(ranLoc$seqnames),
                        start = as.integer(ranLoc$start-1),
                        end = as.integer(ranLoc$end),
                        name = as.integer(1:length(ranLoc$seqnames)),
                        score = rep("NA", length(ranLoc$seqnames)),
                        strand = as.character(ranLoc$strand))
write.table(ranLocbed,
            file = paste0(libName,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, "_randomLoci.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")
