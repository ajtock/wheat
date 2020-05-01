#!/applications/R/R-3.4.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 28.04.2020

# Profile SNP frequency around compartmentalised peaks and random loci

# Wheat subgenome compartments (a.k.a. partitions):
# 1. R1 and R3 (euchromatin)
# 2. R2a and R2b (interstitial)
# 3. C (proximal)
# 4. heterochromatin (interstitial and proximal)
# 5. centromeres (defined by IWGSC (2018) Science 361 using CENH3 ChIP-seq data from Guo et al. (2016) PLOS Genet. 12)

# Usage via Condor submission system on node7:
# csmit -m 50G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./1000exomes_SNP_profiles_around_peaks_commandArgs.R DMC1_Rep1_ChIP euchromatin A 400 2000 2kb 20"

#libName <- "DMC1_Rep1_ChIP"
#region <- "euchromatin"
#genomeName <- "A"
#bodyLength <- 400
#flankSize <- 2000
#flankName <- "2kb"
#winSize <- 20

args <- commandArgs(trailingOnly = T)
libName <- args[1]
region <- args[2]
genomeName <- args[3]
bodyLength <- as.numeric(args[4])
flankSize <- as.numeric(args[5])
flankName <- as.character(args[6])
winSize <- as.numeric(args[7])

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)

matDir <- paste0("matrices/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

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
if(region == "euchromatin") {
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
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
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
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Load peaks in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           libName,
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           genomeName, "genome_", region, ".bed"),
                    header = F)
colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand")
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start+1,
                                    end = peaks$end),
                   strand = peaks$strand,
                   number = peaks$name)
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
ranLoc <- read.table(paste0("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                            libName,
                            "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                            genomeName, "genome_", region, "_randomLoci.bed"),
                     header = F)
colnames(ranLoc) <- c("chr", "start", "end", "name", "score", "strand")
ranLocGR <- GRanges(seqnames = ranLoc$chr,
                    ranges = IRanges(start = ranLoc$start+1,
                                     end = ranLoc$end),
                    strand = ranLoc$strand,
                    number = ranLoc$name)
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

# Load SNPs
SNPs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/He_Akhunov_2019_NatGenet_1000exomes_SNPs/",
                          "all.GP08_mm75_het3_publication01142019.vcf"),
                   header = F, skip = 31,
                   colClasses = c(rep(NA, 2),
                                  "NULL",
                                  rep(NA, 2),
                                  rep("NULL", 2),
                                  NA,
                                  rep("NULL", 812)))
colnames(SNPs) <- c("chr", "pos", "ref", "alt", "info")
## all exome SNPs
SNPsGR <- GRanges(seqnames = SNPs$chr,
                  ranges = IRanges(start = SNPs$pos,
                                   end = SNPs$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs)[1]))

# Define matrix and column mean frequency outfile (mean profiles)
outDF <- list(paste0(matDir,
                     "exome_all_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                     "_matrix_bin", winSize, "bp_flank", flankName, ".tab"),
              paste0(matDir,
                     "exome_all_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                     "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, ".tab"))
outDFcolMeans <- list(paste0(matDir,
                             "exome_all_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                             "_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"),
                      paste0(matDir,
                             "exome_all_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                             "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing frequency values around target and random loci
covMatrix(signal = SNPsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = bodyLength,
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libName, " peaks exome all SNP frequency profile calculation complete"))


## transition
SNPs_transition <- SNPs[(SNPs$ref == "A" | SNPs$ref == "G") & (SNPs$alt == "G" | SNPs$alt == "A") |
                        (SNPs$ref == "C" | SNPs$ref == "T") & (SNPs$alt == "T" | SNPs$alt == "C"),]
SNPsGR <- GRanges(seqnames = SNPs_transition$chr,
                  ranges = IRanges(start = SNPs_transition$pos,
                                   end = SNPs_transition$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs_transition)[1]))

# Define matrix and column mean frequency outfile (mean profiles)
outDF <- list(paste0(matDir,
                     "exome_transition_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                     "_matrix_bin", winSize, "bp_flank", flankName, ".tab"),
              paste0(matDir,
                     "exome_transition_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                     "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, ".tab"))
outDFcolMeans <- list(paste0(matDir,
                             "exome_transition_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                             "_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"),
                      paste0(matDir,
                             "exome_transition_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                             "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing frequency values around target and random loci
covMatrix(signal = SNPsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = bodyLength,
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libName, " peaks exome transition SNP frequency profile calculation complete"))


## transversion
SNPs_transversion <- SNPs[(SNPs$ref == "A" | SNPs$ref == "G") & (SNPs$alt == "C" | SNPs$alt == "T") |
                          (SNPs$ref == "C" | SNPs$ref == "T") & (SNPs$alt == "A" | SNPs$alt == "G"),]
stopifnot((dim(SNPs_transition)[1] +
          dim(SNPs_transversion)[1]) ==
          dim(SNPs)[1])
SNPsGR <- GRanges(seqnames = SNPs_transversion$chr,
                  ranges = IRanges(start = SNPs_transversion$pos,
                                   end = SNPs_transversion$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs_transversion)[1]))

# Define matrix and column mean frequency outfile (mean profiles)
outDF <- list(paste0(matDir,
                     "exome_transversion_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                     "_matrix_bin", winSize, "bp_flank", flankName, ".tab"),
              paste0(matDir,
                     "exome_transversion_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                     "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, ".tab"))
outDFcolMeans <- list(paste0(matDir,
                             "exome_transversion_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                             "_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"),
                      paste0(matDir,
                             "exome_transversion_SNPs_around_DMC1_peaks_in_", genomeName, "genome_", region,
                             "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing frequency values around target and random loci
covMatrix(signal = SNPsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = bodyLength,
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libName, " peaks exome transversion SNP frequency profile calculation complete"))
