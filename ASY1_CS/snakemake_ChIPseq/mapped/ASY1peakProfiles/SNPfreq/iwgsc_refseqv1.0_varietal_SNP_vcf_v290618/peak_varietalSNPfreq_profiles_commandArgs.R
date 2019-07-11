#!/applications/R/R-3.4.0/bin/Rscript

# Profile mean coverage around peaks and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./peak_varietalSNPfreq_profiles_commandArgs.R ASY1_CS ASY1_CS_Rep1_ChIP p0.001_q0.01 2000 2kb 20 'A' 'euchromatin'"

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
#region <- "euchromatin"

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
                           "_noMinWidth_in_", genomeName, "genome_", region, ".bed"),
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
ranLoc <- read.table(paste0(peakDir,
                            libNameChIP,
                            "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel,
                            "_noMinWidth_in_", genomeName, "genome_", region, "_randomLoci.bed"),
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
if(length(mask_ranLoc_overlap) > 1) {
  ranLocGR <- ranLocGR[-subjectHits(mask_ranLoc_overlap)]
}

# Load SNPS
SNPs <- read.table("all_filtered_snps_allaccessions_allploidy_snpeff.vcf",
                   header = F, skip = 6,
                   colClasses = c(rep(NA, 2),
                                  "NULL",
                                  rep(NA, 2),
                                  rep("NULL", 2),
                                  NA,
                                  rep("NULL", 61)))
colnames(SNPs) <- c("chr", "pos", "ref", "alt", "info")
### all SNPs
#SNPsGR <- GRanges(seqnames = SNPs$chr,
#                  ranges = IRanges(start = SNPs$pos,
#                                   end = SNPs$pos),
#                  strand = "*",
#                  coverage = rep(1, dim(SNPs)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_all_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_all_SNP_frequency_ranLoc_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_all_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_all_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#          feature = peaksGR,
#          ranLoc = ranLocGR,
#          featureSize = mean(width(peaksGR)),
#          flankSize = flankSize,
#          winSize = winSize,
#          outDF = outDF,
#          outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks SNP frequency profile calculation complete"))
#
#
### upstream_gene_variant
#SNPs_upstream_gene_variant <- SNPs[grep("upstream_gene_variant", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_upstream_gene_variant$chr,
#                  ranges = IRanges(start = SNPs_upstream_gene_variant$pos,
#                                   end = SNPs_upstream_gene_variant$pos),
#                  strand = "*",
#                  coverage = rep(1, dim(SNPs_upstream_gene_variant)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_upstream_gene_variant_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_upstream_gene_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_upstream_gene_variant_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_upstream_gene_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#          feature = peaksGR,
#          ranLoc = ranLocGR,
#          featureSize = mean(width(peaksGR)),
#          flankSize = flankSize,
#          winSize = winSize,
#          outDF = outDF,
#          outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks upstream_gene_variant SNP frequency profile calculation complete"))
#
#
### downstream_gene_variant
#SNPs_downstream_gene_variant <- SNPs[grep("downstream_gene_variant", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_downstream_gene_variant$chr,
#                  ranges = IRanges(start = SNPs_downstream_gene_variant$pos,
#                                   end = SNPs_downstream_gene_variant$pos),
#                  strand = "*",
#                  coverage = rep(1, dim(SNPs_downstream_gene_variant)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_downstream_gene_variant_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_downstream_gene_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_downstream_gene_variant_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_downstream_gene_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#          feature = peaksGR,
#          ranLoc = ranLocGR,
#          featureSize = mean(width(peaksGR)),
#          flankSize = flankSize,
#          winSize = winSize,
#          outDF = outDF,
#          outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks downstream_gene_variant SNP frequency profile calculation complete"))
#
#
### missense_variant
#SNPs_missense_variant <- SNPs[grep("missense_variant", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_missense_variant$chr,
#                  ranges = IRanges(start = SNPs_missense_variant$pos,
#                                   end = SNPs_missense_variant$pos),
#                  strand = "*",
#                  coverage = rep(1, dim(SNPs_missense_variant)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_missense_variant_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_missense_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_missense_variant_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_missense_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#          feature = peaksGR,
#          ranLoc = ranLocGR,
#          featureSize = mean(width(peaksGR)),
#          flankSize = flankSize,
#          winSize = winSize,
#          outDF = outDF,
#          outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks missense_variant SNP frequency profile calculation complete"))
#
#
### synonymous_variant
#SNPs_synonymous_variant <- SNPs[grep("synonymous_variant", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_synonymous_variant$chr,
#                  ranges = IRanges(start = SNPs_synonymous_variant$pos,
#                                   end = SNPs_synonymous_variant$pos),
#                  strand = "*",
#                  coverage = rep(1, dim(SNPs_synonymous_variant)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_synonymous_variant_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_synonymous_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_synonymous_variant_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_synonymous_variant_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#          feature = peaksGR,
#          ranLoc = ranLocGR,
#          featureSize = mean(width(peaksGR)),
#          flankSize = flankSize,
#          winSize = winSize,
#          outDF = outDF,
#          outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks synonymous_variant SNP frequency profile calculation complete"))
#
#
### HIGH
#SNPs_HIGH <- SNPs[grep("HIGH", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_HIGH$chr,
#                  ranges = IRanges(start = SNPs_HIGH$pos,
#                                   end = SNPs_HIGH$pos),
#                  strand = "*",
#                  coverage = rep(1, dim(SNPs_HIGH)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_HIGH_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#	      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_HIGH_SNP_frequency_ranLoc_smoothed_target_and_",
#                      flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_HIGH_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#		      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_HIGH_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#          feature = peaksGR,
#          ranLoc = ranLocGR,
#          featureSize = mean(width(peaksGR)),
#          flankSize = flankSize,
#          winSize = winSize,
#          outDF = outDF,
#          outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks HIGH SNP frequency profile calculation complete"))
#
#
### MODERATE
#SNPs_MODERATE <- SNPs[grep("MODERATE", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_MODERATE$chr,
#		  ranges = IRanges(start = SNPs_MODERATE$pos,
#				   end = SNPs_MODERATE$pos),
#		  strand = "*",
#		  coverage = rep(1, dim(SNPs_MODERATE)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_MODERATE_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#	      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_MODERATE_SNP_frequency_ranLoc_smoothed_target_and_",
#                      flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_MODERATE_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#		      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_MODERATE_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#	  feature = peaksGR,
#	  ranLoc = ranLocGR,
#	  featureSize = mean(width(peaksGR)),
#	  flankSize = flankSize,
#	  winSize = winSize,
#	  outDF = outDF,
#	  outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks MODERATE SNP frequency profile calculation complete"))
#
#
### LOW
#SNPs_LOW <- SNPs[grep("LOW", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_LOW$chr,
#		  ranges = IRanges(start = SNPs_LOW$pos,
#				   end = SNPs_LOW$pos),
#		  strand = "*",
#		  coverage = rep(1, dim(SNPs_LOW)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_LOW_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#	      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_LOW_SNP_frequency_ranLoc_smoothed_target_and_",
#                      flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_LOW_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#		      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_LOW_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#	  feature = peaksGR,
#	  ranLoc = ranLocGR,
#	  featureSize = mean(width(peaksGR)),
#	  flankSize = flankSize,
#	  winSize = winSize,
#	  outDF = outDF,
#	  outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks LOW SNP frequency profile calculation complete"))
#
#
### MODIFIER
#SNPs_MODIFIER <- SNPs[grep("MODIFIER", SNPs$info),]
#SNPsGR <- GRanges(seqnames = SNPs_MODIFIER$chr,
#		  ranges = IRanges(start = SNPs_MODIFIER$pos,
#				   end = SNPs_MODIFIER$pos),
#		  strand = "*",
#		  coverage = rep(1, dim(SNPs_MODIFIER)[1]))
#
## Define matrix and column mean coverage outfile (mean profiles)
#outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_MODIFIER_SNP_frequency_feature_smoothed_target_and_",
#                     flankName, "_flank_dataframe.txt"),
#	      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                     "_MODIFIER_SNP_frequency_ranLoc_smoothed_target_and_",
#                      flankName, "_flank_dataframe.txt"))
#outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_MODIFIER_SNP_frequency_feature_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"),
#		      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
#                             "_MODIFIER_SNP_frequency_ranLoc_smoothed_target_and_",
#                             flankName, "_flank_dataframe_colMeans.txt"))
#
## Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
#covMatrix(signal = SNPsGR,
#	  feature = peaksGR,
#	  ranLoc = ranLocGR,
#	  featureSize = mean(width(peaksGR)),
#	  flankSize = flankSize,
#	  winSize = winSize,
#	  outDF = outDF,
#	  outDFcolMeans = outDFcolMeans)
#print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks MODIFIER SNP frequency profile calculation complete"))


## transition
SNPs_transition <- SNPs[(SNPs$ref == "A" | SNPs$ref == "G") & (SNPs$alt == "G" | SNPs$alt == "A") |
                        (SNPs$ref == "C" | SNPs$ref == "T") & (SNPs$alt == "T" | SNPs$alt == "C"),]
SNPsGR <- GRanges(seqnames = SNPs_transition$chr,
		  ranges = IRanges(start = SNPs_transition$pos,
				   end = SNPs_transition$pos),
		  strand = "*",
		  coverage = rep(1, dim(SNPs_transition)[1]))

# Define matrix and column mean coverage outfile (mean profiles)
outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_transition_SNP_frequency_feature_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"),
	      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_transition_SNP_frequency_ranLoc_smoothed_target_and_",
                      flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_transition_SNP_frequency_feature_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"),
		      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_transition_SNP_frequency_ranLoc_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing normalised coverage values around target and random loci
covMatrix(signal = SNPsGR,
	  feature = peaksGR,
	  ranLoc = ranLocGR,
	  featureSize = mean(width(peaksGR)),
	  flankSize = flankSize,
	  winSize = winSize,
	  outDF = outDF,
	  outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks transition SNP frequency profile calculation complete"))


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

# Define matrix and column mean coverage outfile (mean profiles)
outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_transversion_SNP_frequency_feature_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"),
              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_transversion_SNP_frequency_ranLoc_smoothed_target_and_",
                      flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_transversion_SNP_frequency_feature_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_transversion_SNP_frequency_ranLoc_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing normalised coverage values around target and random loci
covMatrix(signal = SNPsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(peaksGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks transversion SNP frequency profile calculation complete"))


## intron_variant
SNPs_intron_variant <- SNPs[grep("intron_variant", SNPs$info),]
SNPsGR <- GRanges(seqnames = SNPs_intron_variant$chr,
                  ranges = IRanges(start = SNPs_intron_variant$pos,
                                   end = SNPs_intron_variant$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs_intron_variant)[1]))

# Define matrix and column mean coverage outfile (mean profiles)
outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_intron_variant_SNP_frequency_feature_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"),
              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_intron_variant_SNP_frequency_ranLoc_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_intron_variant_SNP_frequency_feature_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_intron_variant_SNP_frequency_ranLoc_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing normalised coverage values around target and random loci
covMatrix(signal = SNPsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(peaksGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks intron_variant SNP frequency profile calculation complete"))


## intergenic
SNPs_intergenic <- SNPs[grep("intergenic", SNPs$info),]
SNPsGR <- GRanges(seqnames = SNPs_intergenic$chr,
                  ranges = IRanges(start = SNPs_intergenic$pos,
                                   end = SNPs_intergenic$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs_intergenic)[1]))

# Define matrix and column mean coverage outfile (mean profiles)
outDF <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_intergenic_SNP_frequency_feature_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"),
              paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                     "_intergenic_SNP_frequency_ranLoc_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_intergenic_SNP_frequency_feature_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, libNameChIP, "_peaks_in_", genomeName, "genome_", region,
                             "_intergenic_SNP_frequency_ranLoc_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing normalised coverage values around target and random loci
covMatrix(signal = SNPsGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(peaksGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", libNameChIP, " peaks intergenic SNP frequency profile calculation complete"))
