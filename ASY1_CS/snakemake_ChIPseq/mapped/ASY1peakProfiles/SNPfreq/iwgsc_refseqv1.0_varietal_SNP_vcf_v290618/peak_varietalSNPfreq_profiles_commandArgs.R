#!/applications/R/R-3.4.0/bin/Rscript

# Profile mean coverage around peaks and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "./peak_varietalSNPfreq_profiles_commandArgs.R 2000 2kb 20 'A' 'euchromatin'"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)

#flankSize <- 2000
#flankName <- "2kb"
#winSize <- 20
#genomeName <- "A"
#region <- "euchromatin"

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
genomeName <- args[4]
region <- args[5]

matDir <- paste0("matrices/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
} else if(region == "pericentromeres") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
} else if(region == "pericentromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
} else if(region == "genomewide") {
  maskGR <- GRanges()
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Load peaks in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           "ASY1_CS_Rep1_ChIP",
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           genomeName, "genome_", region, ".bed"),
                    header = F)
colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand")
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start+1,
                                    end = peaks$end),
                   strand = peaks$strand,
                   peakID = peaks$name)
peaksGR <- peaksGR[seqnames(peaksGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
if(length(mask_peaks_overlap) > 0) {
  peaksGR <- peaksGR[-subjectHits(mask_peaks_overlap)]
}
# Load ranLoc in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLoc <- read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                            "ASY1_CS_Rep1_ChIP",
                            "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                            genomeName, "genome_", region, "_randomLoci.bed"),
                     header = F)
colnames(ranLoc) <- c("chr", "start", "end", "name", "score", "strand")
ranLocGR <- GRanges(seqnames = ranLoc$chr,
                   ranges = IRanges(start = ranLoc$start+1,
                                    end = ranLoc$end),
                   strand = ranLoc$strand,
                   peakID = ranLoc$name)
ranLocGR <- ranLocGR[seqnames(ranLocGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_ranLoc_overlap <- findOverlaps(query = maskGR,
                                   subject = ranLocGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
if(length(mask_ranLoc_overlap) > 0) {
  ranLocGR <- ranLocGR[-subjectHits(mask_ranLoc_overlap)]
}

# Load SNPS
SNPs <- read.table("all_filtered_snps_allaccessions_allploidy_snpeff.vcf",
                   header = F, skip = 6,
                   colClasses = c(rep(NA, 2), rep("NULL", 67)))
colnames(SNPs) <- c("chr", "pos")
SNPsGR <- GRanges(seqnames = SNPs$chr,
                  ranges = IRanges(start = SNPs$pos,
                                   end = SNPs$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs)[1]))

# Define matrix and column mean coverage outfile (mean profiles)
outDF <- list(paste0(matDir, "ASY1_CS_peaks_in_", genomeName, "genome_", region,
                     "_SNP_frequency_feature_smoothed_target_and_", flankName, "_flank_dataframe.txt"),
              paste0(matDir, "ASY1_CS_peaks_in_", genomeName, "genome_", region,
                     "_SNP_frequency_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, "ASY1_CS_peaks_in_", genomeName, "genome_", region,
                             "_SNP_frequency_feature_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, "ASY1_CS_peaks_in_", genomeName, "genome_", region,
                             "_SNP_frequency_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
covMatrix(signal = SNPsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(peaksGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ASY1 peaks SNP frequency profile calculation complete"))
