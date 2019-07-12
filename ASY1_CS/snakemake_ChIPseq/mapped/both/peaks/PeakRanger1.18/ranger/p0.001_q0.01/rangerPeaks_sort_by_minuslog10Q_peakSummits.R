#!/applications/R/R-3.3.2/bin/Rscript

# Define peak coordinates as peak summit-100 bp to peak summit+100 bp
# so that all peaks have a common width of 201 bp for motif analysis
# Generate peak gff and bed files (0-based start coordinates)

# Usage:
# ./rangerPeaks_sort_by_minuslog10Q_peakSummits.R ASY1_CS_Rep1_ChIP 'euchromatin' 'A'

source("/projects/ajt200/Rfunctions/locus_midpoint.R")
library(GenomicRanges)

args <- commandArgs(trailingOnly = T)
libName <- args[1]
region <- args[2]
genomeName <- args[3]

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
peaks <- read.table(paste0(libName,
                           "_peaks_peakranger_ranger_p0.001_q0.01_log2TreadsNormCreads.narrowPeak"))
peaks <- cbind(peaks[,1:3],
               peaks[,7:10])
colnames(peaks) <- c("chr", "start0based", "end",
                     "sigval", "pval", "qval", "summit0based")
peaks <- data.frame(chr = peaks$chr,
                    start = peaks$start0based+1,
                    end = peaks$end,
                    sigval = peaks$sigval,
                    pval = peaks$pval,
                    qval = peaks$qval,
                    summit0based = peaks$summit0based)
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start,
                                    end = peaks$end),
                   strand = "*",
                   sigval = peaks$sigval,
                   pval = peaks$pval,
                   qval = peaks$qval,
                   summit0based = peaks$summit0based)
peaksGR <- peaksGR[grep(genomeName,
                        seqnames(peaksGR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaksGR <- peaksGR[-subjectHits(mask_peaks_overlap)]
print("***********peaks***********")
print(peaksGR)

# Sort by decreasing -log10(qval)
peaksGR <- sort(peaksGR, by = ~ qval, decreasing = T)


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
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(as.vector(ranges(regionChrGR))),
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
                       feature = as.character(rep(paste0(libName, "_peak"))),
                       start = as.integer(peaks$start),
                       end = as.integer(peaks$end),
                       qval = as.numeric(peaks$qval),
                       strand = as.character(rep(".")),
                       frame = as.character(rep(".")),
                       summit0based = as.integer(peaks$summit0based))
write.table(peaksgff,
            file = paste0(libName,
                          "_rangerPeaksGR_minuslog10Qsorted_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, ".gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")

peaksbed <- data.frame(chr = as.character(peaks$seqnames),
                       start = as.integer(peaks$start-1),
                       end = as.integer(peaks$end),
                       name = as.integer(1:length(peaks$seqnames)),
                       qval = as.numeric(peaks$qval),
                       summit0based = as.integer(peaks$summit0based))
write.table(peaksbed,
            file = paste0(libName,
                          "_rangerPeaksGR_minuslog10Qsorted_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, ".bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# Extract peak summits +/- 200 bp for use in weeder2 motif analysis
peaksGR_summits <- GRanges(seqnames = seqnames(peaksGR),
                           ranges = IRanges(start = start(peaksGR)+peaksGR$summit0based,
                                            end = start(peaksGR)+peaksGR$summit0based),
                           strand = "*",
                           sigval = peaksGR$sigval,
                           pval = peaksGR$pval,
                           qval = peaksGR$qval,
                           summit0based = peaksGR$summit0based)
peaksGR_summits200bp <- locMidpointFlank(x = peaksGR_summits,
                                         leftFlank = 100,
                                         rightFlank = 100)
peaks_summits200bp_bed <- data.frame(chr = seqnames(peaksGR_summits200bp),
                                     start = start(peaksGR_summits200bp)-1,
                                     end = end(peaksGR_summits200bp))
write.table(peaks_summits200bp_bed,
            file = paste0(region, "/motifs_summits200bp/", libName,
                          "_rangerPeaksGR_minuslog10Qsorted_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, "_summits200bp.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# ranLoc
ranLoc <- data.frame(ranLocGR)

ranLocgff <- data.frame(chr = as.character(ranLoc$seqnames),
                        source = as.character(rep(".")),
                        feature = as.character(rep(paste0(libName, "_peak_ranLoc"))),
                        start = as.integer(ranLoc$start),
                        end = as.integer(ranLoc$end),
                        score = as.character(rep(".")),
                        strand = as.character(rep(".")),
                        frame = as.character(rep(".")),
                        attribute = as.character(rep(".")))
write.table(ranLocgff,
            file = paste0(libName,
                          "_rangerPeaksGR_minuslog10Qsorted_p0.001_q0.01_noMinWidth_in_",
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
                          "_rangerPeaksGR_minuslog10Qsorted_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, "_randomLoci.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# Extract ranLoc summits +/- 200 bp for use in weeder2 motif analysis
ranLocGR_summits <- GRanges(seqnames = seqnames(ranLocGR),
                            ranges = IRanges(start = start(ranLocGR)+round((end(ranLocGR)-start(ranLocGR))/2),
                                             end = start(ranLocGR)+round((end(ranLocGR)-start(ranLocGR))/2)),
                            strand = "*")
ranLocGR_summits200bp <- locMidpointFlank(x = ranLocGR_summits,
                                          leftFlank = 100,
                                          rightFlank = 100)
ranLoc_summits200bp_bed <- data.frame(chr = seqnames(ranLocGR_summits200bp),
                                      start = start(ranLocGR_summits200bp)-1,
                                      end = end(ranLocGR_summits200bp))
write.table(ranLoc_summits200bp_bed,
            file = paste0(region, "/motifs_summits200bp/", libName,
                          "_rangerPeaksGR_minuslog10Qsorted_p0.001_q0.01_noMinWidth_in_",
                          genomeName, "genome_", region, "_summits200bp_randomLoci.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

