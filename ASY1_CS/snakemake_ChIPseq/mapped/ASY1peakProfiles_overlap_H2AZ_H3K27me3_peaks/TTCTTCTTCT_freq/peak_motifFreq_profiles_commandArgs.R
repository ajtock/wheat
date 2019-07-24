#!/applications/R/R-3.4.0/bin/Rscript

# Profile motif frequency around peaks and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./peak_motifFreq_profiles_commandArgs.R ASY1_CS ASY1_CS_Rep1_ChIP p0.001_q0.01 2000 2kb 20 'A' 'heterochromatin'"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)

#markChIP <- "ASY1_CS"
#libNameChIP <- "ASY1_CS_Rep1_ChIP"
#sigLevel <- "p0.001_q0.01"
#flankSize <- 2000
#flankName <- "2kb"
#winSize <- 20
#genomeName <- "A"
#region <- "heterochromatin"

args <- commandArgs(trailingOnly = T)
markChIP <- args[1]
libNameChIP <- args[2]
sigLevel <- args[3]
flankSize <- as.numeric(args[4])
flankName <- as.character(args[5])
winSize <- as.numeric(args[6])
genomeName <- args[7]
region <- args[8]

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
if(libNameChIP %in% c("H3K4me3_ChIP_SRR6350668",
                      "H3K27me3_ChIP_SRR6350666",
                      "H3K36me3_ChIP_SRR6350670",
                      "H3K9ac_ChIP_SRR6350667",
                      "CENH3_ChIP_SRR1686799")) {
  peakDir <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                    markChIP,
                    "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                    sigLevel, "/")
} else {
  peakDir <- paste0("/home/ajt200/analysis/wheat/",
                    markChIP,
                    "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                    sigLevel, "/")
}
peaks <- read.table(paste0(peakDir,
                           libNameChIP,
                           "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel,
                           "_noMinWidth_in_", genomeName, "genome_", region,
                           "_overlapping_H2AZ_Rep1_ChIP_and_H3K27me3_ChIP_SRR6350666_peaks.bed"),
                    header = F)
colnames(peaks) <- c("chr", "start", "end", "name", "qval", "summit0based")
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start+1,
                                    end = peaks$end),
                   strand = "*",
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
ranLoc <- read.table(paste0(peakDir,
                            libNameChIP,
                            "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel,
                            "_noMinWidth_in_", genomeName, "genome_", region,
                            "_overlapping_H2AZ_Rep1_ChIP_and_H3K27me3_ChIP_SRR6350666_peaks_randomLoci.bed"),
                     header = F)
colnames(ranLoc) <- c("chr", "start", "end", "name", "qval", "summit0based")
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
if(length(mask_ranLoc_overlap) > 1) {
  ranLocGR <- ranLocGR[-subjectHits(mask_ranLoc_overlap)]
}

# Load motif-matching loci
load(paste0("/home/ajt200/analysis/wheat/recombination_rate/ASY1peakOverlap_H2AZ_H3K27me3_CTTrepeats/CTTrepeats/",
            region, "/",
            "MAT1_TTCTTCTTCT_motif_TTCTTCTTCT_matchPWM_loci_in_",
            genomeName, "genome_", region, ".RData"))
strand(hitsGR) <- "*"
hitsGR$coverage <- rep(1, length(hitsGR))

# Define matrix and column mean outfiles (mean profiles)
outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_motif_frequency_feature_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"),
              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_motif_frequency_ranLoc_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_motif_frequency_feature_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_motif_frequency_ranLoc_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing normalised coverage values around target and random loci
covMatrix(signal = hitsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(peaksGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks motif frequency profile calculation complete"))
