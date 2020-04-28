#!/applications/R/R-3.3.2/bin/Rscript

# Obtain ASY1 peaks that overlap H2A.Z and H3K27me3 peaks,
# to be ordered by decreasing ASY1 coverage and analysed for motif enrichment
# Generate peak gff and bed files (0-based start coordinates)

# Usage:
# ./rangerPeaks_ASY1_H2AZ_H3K27me3_peakOverlap.R ASY1_CS ASY1_CS_Rep1_ChIP p0.001_q0.01 H2AZ H2AZ_Rep1_ChIP p0.05_q0.05 H3K27me3 H3K27me3_ChIP_SRR6350666 p0.05_q0.05 'euchromatin' 'A'

library(GenomicRanges)

#mark1 <- "ASY1_CS"
#lib1 <- "ASY1_CS_Rep1_ChIP"
#sigLevel1 <- "p0.001_q0.01"
#mark2 <- "H2AZ"
#lib2 <- "H2AZ_Rep1_ChIP"
#sigLevel2 <- "p0.05_q0.05"
#mark3 <- "H3K27me3"
#lib3 <- "H3K27me3_ChIP_SRR6350666"
#sigLevel3 <- "p0.05_q0.05"
#region <- "euchromatin"
#genomeName <- "A"

args <- commandArgs(trailingOnly = T)
mark1 <- args[1]
lib1 <- args[2]
sigLevel1 <- args[3]
mark2 <- args[4]
lib2 <- args[5]
sigLevel2 <- args[6]
mark3 <- args[7]
lib3 <- args[8]
sigLevel3 <- args[9]
region <- args[10]
genomeName <- args[11]

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
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Import peaks
# lib1
if(lib1 %in% c("H3K4me3_ChIP_SRR6350668",
               "H3K27me3_ChIP_SRR6350666",
               "H3K36me3_ChIP_SRR6350670",
               "H3K9ac_ChIP_SRR6350667",
               "CENH3_ChIP_SRR1686799")) {
  dir1 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/", mark1,
                 "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                 sigLevel1, "/")
} else {
  dir1 <- paste0("/home/ajt200/analysis/wheat/", mark1,
                 "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                 sigLevel1, "/")
}
load(paste0(dir1, lib1, "_rangerPeaksGRmergedOverlaps_minuslog10_",
            sigLevel1, "_noMinWidth.RData"))
peaks1GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
peaks1GR <- peaks1GR[grep(genomeName,
                          seqnames(peaks1GR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaks1GR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaks1GR <- peaks1GR[-subjectHits(mask_peaks_overlap)]
print("***********peaks1***********")
print(peaks1GR)

# Append FDRs (qvals) and peak summits as columns
peaks1_qvalSorted <- read.table(paste0(dir1, lib1,
                                       "_rangerPeaksGR_minuslog10Qsorted_",
                                       sigLevel1, "_noMinWidth_in_",
                                       genomeName, "genome_", region, ".gff"))
chrs <- chrs[grep(genomeName, chrs)]
peaks1_qval <- data.frame()
for(i in seq_along(chrs)) {
  peaks1_qval_chr <- peaks1_qvalSorted[peaks1_qvalSorted$V1 == chrs[i],]
  peaks1GR_chr <- peaks1GR[seqnames(peaks1GR) == chrs[i]]
  peaks1_qval_chr <- peaks1_qval_chr[peaks1_qval_chr$V4 %in% start(peaks1GR_chr),]
  peaks1_qval_chrSorted <- peaks1_qval_chr[order(peaks1_qval_chr$V4),]
  peaks1_qval <- rbind(peaks1_qval, peaks1_qval_chrSorted)
}

peaks1GR <- GRanges(peaks1GR,
                    qval = peaks1_qval$V6,
                    summit0based = peaks1_qval$V9) 

# lib2
if(lib2 %in% c("H3K4me3_ChIP_SRR6350668",
               "H3K27me3_ChIP_SRR6350666",
               "H3K36me3_ChIP_SRR6350670",
               "H3K9ac_ChIP_SRR6350667",
               "CENH3_ChIP_SRR1686799")) {
  dir2 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/", mark2,
                 "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                 sigLevel2, "/")
} else {
  dir2 <- paste0("/home/ajt200/analysis/wheat/", mark2,
                 "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                 sigLevel2, "/")
}
load(paste0(dir2, lib2, "_rangerPeaksGRmergedOverlaps_minuslog10_",
            sigLevel2, "_noMinWidth.RData"))
peaks2GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
peaks2GR <- peaks2GR[grep(genomeName,
                          seqnames(peaks2GR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaks2GR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaks2GR <- peaks2GR[-subjectHits(mask_peaks_overlap)]
print("***********peaks2***********")
print(peaks2GR)

# lib3
if(lib3 %in% c("H3K4me3_ChIP_SRR6350668",
               "H3K27me3_ChIP_SRR6350666",
               "H3K36me3_ChIP_SRR6350670",
               "H3K9ac_ChIP_SRR6350667",
               "CENH3_ChIP_SRR1686799")) {
  dir3 <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/", mark3,
                 "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                 sigLevel3, "/")
} else {
  dir3 <- paste0("/home/ajt200/analysis/wheat/", mark3,
                 "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                 sigLevel3, "/")
}
load(paste0(dir3, lib3, "_rangerPeaksGRmergedOverlaps_minuslog10_",
            sigLevel3, "_noMinWidth.RData"))
peaks3GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
peaks3GR <- peaks3GR[grep(genomeName,
                          seqnames(peaks3GR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaks3GR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaks3GR <- peaks3GR[-subjectHits(mask_peaks_overlap)]
print("***********peaks3***********")
print(peaks3GR)

# Obtain overlapping peaks
peaks1_peaks2_overlap <- findOverlaps(query = peaks1GR,
                                      subject = peaks2GR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
peaks1_peaks2_overlapGR <- unique(peaks1GR[queryHits(peaks1_peaks2_overlap)])

peaks1_peaks3_overlap <- findOverlaps(query = peaks1GR,
                                      subject = peaks3GR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
peaks1_peaks3_overlapGR <- unique(peaks1GR[queryHits(peaks1_peaks3_overlap)])

peaks1_peaks2_peaks3_overlap <- findOverlaps(query = peaks1_peaks2_overlapGR,
                                             subject = peaks3GR,
                                             type = "any",
                                             select = "all",
                                             ignore.strand = TRUE)
peaks1_peaks2_peaks3_overlapGR <- unique(peaks1_peaks2_overlapGR[queryHits(peaks1_peaks2_peaks3_overlap)])


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
ranLoc1GR <- GRanges()
for(i in 1:length(chrs)) {
  peaksChrGR <- peaks1_peaks2_overlapGR[seqnames(peaks1_peaks2_overlapGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(peaksChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(as.vector(ranges(regionChrGR))),
                                      n = length(peaksChrGR))
  ranLocChrIR <- IRanges(start = ranLocChrStart,
                         width = width(peaksChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = ranLocChrIR,
                         strand = "*")
  ranLoc1GR <- append(ranLoc1GR, ranLocChrGR)
}

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as peaksGR
chrs <- chrs[grep(genomeName, chrs)]
ranLoc2GR <- GRanges()
for(i in 1:length(chrs)) {
  peaksChrGR <- peaks1_peaks3_overlapGR[seqnames(peaks1_peaks3_overlapGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(peaksChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(as.vector(ranges(regionChrGR))),
                                      n = length(peaksChrGR))
  ranLocChrIR <- IRanges(start = ranLocChrStart,
                         width = width(peaksChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = ranLocChrIR,
                         strand = "*")
  ranLoc2GR <- append(ranLoc2GR, ranLocChrGR)
}

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as peaksGR
chrs <- chrs[grep(genomeName, chrs)]
ranLoc3GR <- GRanges()
for(i in 1:length(chrs)) {
  peaksChrGR <- peaks1_peaks2_peaks3_overlapGR[seqnames(peaks1_peaks2_peaks3_overlapGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(peaksChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(as.vector(ranges(regionChrGR))),
                                      n = length(peaksChrGR))
  ranLocChrIR <- IRanges(start = ranLocChrStart,
                         width = width(peaksChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = ranLocChrIR,
                         strand = "*")
  ranLoc3GR <- append(ranLoc3GR, ranLocChrGR)
}

# Convert into GFF3 and BED formats
# peaks1_peaks2_OL
peaks1_peaks2_OL <- data.frame(peaks1_peaks2_overlapGR)

peaks1_peaks2_OLgff <- data.frame(chr = as.character(peaks1_peaks2_OL$seqnames),
                                  source = as.character(rep(".")),
                                  feature = as.character(rep(paste0(lib1, "_peak"))),
                                  start = as.integer(peaks1_peaks2_OL$start),
                                  end = as.integer(peaks1_peaks2_OL$end),
                                  qval = as.numeric(peaks1_peaks2_OL$qval),
                                  strand = as.character(rep(".")),
                                  phase = as.character(rep(".")),
                                  summit0based = as.integer(peaks1_peaks2_OL$summit0based))
write.table(peaks1_peaks2_OLgff,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_peaks.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

peaks1_peaks2_OLbed <- data.frame(chr = as.character(peaks1_peaks2_OL$seqnames),
                                  start = as.integer(peaks1_peaks2_OL$start-1),
                                  end = as.integer(peaks1_peaks2_OL$end),
                                  name = as.integer(1:length(peaks1_peaks2_OL$seqnames)),
                                  qval = as.numeric(peaks1_peaks2_OL$qval),
                                  summit0based = as.integer(peaks1_peaks2_OL$summit0based))
write.table(peaks1_peaks2_OLbed,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_peaks.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# peaks1_peaks3_OL
peaks1_peaks3_OL <- data.frame(peaks1_peaks3_overlapGR)

peaks1_peaks3_OLgff <- data.frame(chr = as.character(peaks1_peaks3_OL$seqnames),
                                  source = as.character(rep(".")),
                                  feature = as.character(rep(paste0(lib1, "_peak"))),
                                  start = as.integer(peaks1_peaks3_OL$start),
                                  end = as.integer(peaks1_peaks3_OL$end),
                                  qval = as.numeric(peaks1_peaks3_OL$qval),
                                  strand = as.character(rep(".")),
                                  frame = as.character(rep(".")),
                                  summit0based = as.integer(peaks1_peaks3_OL$summit0based))
write.table(peaks1_peaks3_OLgff,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib3, "_peaks.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

peaks1_peaks3_OLbed <- data.frame(chr = as.character(peaks1_peaks3_OL$seqnames),
                                  start = as.integer(peaks1_peaks3_OL$start-1),
                                  end = as.integer(peaks1_peaks3_OL$end),
                                  name = as.integer(1:length(peaks1_peaks3_OL$seqnames)),
                                  qval = as.numeric(peaks1_peaks3_OL$qval),
                                  summit0based = as.integer(peaks1_peaks3_OL$summit0based))
write.table(peaks1_peaks3_OLbed,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib3, "_peaks.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# peaks1_peaks2_peaks3_OL
peaks1_peaks2_peaks3_OL <- data.frame(peaks1_peaks2_peaks3_overlapGR)

peaks1_peaks2_peaks3_OLgff <- data.frame(chr = as.character(peaks1_peaks2_peaks3_OL$seqnames),
                                         source = as.character(rep(".")),
                                         feature = as.character(rep(paste0(lib1, "_peak"))),
                                         start = as.integer(peaks1_peaks2_peaks3_OL$start),
                                         end = as.integer(peaks1_peaks2_peaks3_OL$end),
                                         qval = as.numeric(peaks1_peaks2_peaks3_OL$qval),
                                         strand = as.character(rep(".")),
                                         frame = as.character(rep(".")),
                                         summit0based = as.integer(peaks1_peaks2_peaks3_OL$summit0based))
write.table(peaks1_peaks2_peaks3_OLgff,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_and_", lib3, "_peaks.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

peaks1_peaks2_peaks3_OLbed <- data.frame(chr = as.character(peaks1_peaks2_peaks3_OL$seqnames),
                                         start = as.integer(peaks1_peaks2_peaks3_OL$start-1),
                                         end = as.integer(peaks1_peaks2_peaks3_OL$end),
                                         name = as.integer(1:length(peaks1_peaks2_peaks3_OL$seqnames)),
                                         qval = as.numeric(peaks1_peaks2_peaks3_OL$qval),
                                         summit0based = as.integer(peaks1_peaks2_peaks3_OL$summit0based))
write.table(peaks1_peaks2_peaks3_OLbed,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_and_", lib3, "_peaks.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")


# ranLoc1
ranLoc1 <- data.frame(ranLoc1GR)

ranLoc1gff <- data.frame(chr = as.character(ranLoc1$seqnames),
                         source = as.character(rep(".")),
                         feature = as.character(rep(paste0(lib1, "_peak_ranLoc"))),
                         start = as.integer(ranLoc1$start),
                         end = as.integer(ranLoc1$end),
                         score = as.character(rep(".")),
                         strand = as.character(rep(".")),
                         phase = as.character(rep(".")),
                         attribute = as.character(rep(".")))
write.table(ranLoc1gff,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_peaks_randomLoci.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

ranLoc1bed <- data.frame(chr = as.character(ranLoc1$seqnames),
                         start = as.integer(ranLoc1$start-1),
                         end = as.integer(ranLoc1$end),
                         name = as.integer(1:length(ranLoc1$seqnames)),
                         score = rep("NA", length(ranLoc1$seqnames)),
                         strand = as.character(ranLoc1$strand))
write.table(ranLoc1bed,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_peaks_randomLoci.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# ranLoc2
ranLoc2 <- data.frame(ranLoc2GR)

ranLoc2gff <- data.frame(chr = as.character(ranLoc2$seqnames),
                         source = as.character(rep(".")),
                         feature = as.character(rep(paste0(lib1, "_peak_ranLoc"))),
                         start = as.integer(ranLoc2$start),
                         end = as.integer(ranLoc2$end),
                         score = as.character(rep(".")),
                         strand = as.character(rep(".")),
                         phase = as.character(rep(".")),
                         attribute = as.character(rep(".")))
write.table(ranLoc2gff,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib3, "_peaks_randomLoci.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

ranLoc2bed <- data.frame(chr = as.character(ranLoc2$seqnames),
                         start = as.integer(ranLoc2$start-1),
                         end = as.integer(ranLoc2$end),
                         name = as.integer(1:length(ranLoc2$seqnames)),
                         score = rep("NA", length(ranLoc2$seqnames)),
                         strand = as.character(ranLoc2$strand))
write.table(ranLoc2bed,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib3, "_peaks_randomLoci.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# ranLoc3
ranLoc3 <- data.frame(ranLoc3GR)

ranLoc3gff <- data.frame(chr = as.character(ranLoc3$seqnames),
                         source = as.character(rep(".")),
                         feature = as.character(rep(paste0(lib1, "_peak_ranLoc"))),
                         start = as.integer(ranLoc3$start),
                         end = as.integer(ranLoc3$end),
                         score = as.character(rep(".")),
                         strand = as.character(rep(".")),
                         phase = as.character(rep(".")),
                         attribute = as.character(rep(".")))
write.table(ranLoc3gff,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_and_", lib3, "_peaks_randomLoci.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

ranLoc3bed <- data.frame(chr = as.character(ranLoc3$seqnames),
                         start = as.integer(ranLoc3$start-1),
                         end = as.integer(ranLoc3$end),
                         name = as.integer(1:length(ranLoc3$seqnames)),
                         score = rep("NA", length(ranLoc3$seqnames)),
                         strand = as.character(ranLoc3$strand))
write.table(ranLoc3bed,
            file = paste0(dir1, lib1,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_",
                          sigLevel1, "_noMinWidth_in_",
                          genomeName, "genome_", region,
                          "_overlapping_", lib2, "_and_", lib3, "_peaks_randomLoci.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")
