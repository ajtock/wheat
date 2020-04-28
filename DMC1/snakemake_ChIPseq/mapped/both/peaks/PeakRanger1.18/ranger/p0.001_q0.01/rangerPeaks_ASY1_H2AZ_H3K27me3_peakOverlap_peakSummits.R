#!/applications/R/R-3.5.0/bin/Rscript

# Define peak coordinates as peak summit-(size/2) bp to peak summit+(size/2) bp
# so that all peaks have a common width of size+1 bp for motif analysis
# Generate peak gff and bed files (0-based start coordinates)

# Usage:
# ./rangerPeaks_ASY1_H2AZ_H3K27me3_peakOverlap_peakSummits.R ASY1_CS ASY1_CS_Rep1_ChIP ASY1 ASY1_CS_Rep1_ChIP p0.001_q0.01 'euchromatin' 'A' 200

source("/projects/ajt200/Rfunctions/locus_midpoint.R")
library(GenomicRanges)

#markName <- "ASY1_CS"
#libName <- "ASY1_CS_Rep1_ChIP"
#featureName <- "ASY1"
#featureLibName <- "ASY1_CS_Rep1_ChIP"
#sigLevel <- "p0.001_q0.01"
#region <- "euchromatin"
#genomeName <- "A"
#size <- 200

args <- commandArgs(trailingOnly = T)
markName <- args[1]
libName <- args[2]
featureName <- args[3]
featureLibName <- args[4]
sigLevel <- args[5]
region <- args[6]
genomeName <- args[7]
size <- as.numeric(args[8])

regionDir <- paste0(region, "/")
motifDir <- paste0(regionDir,
                   "motifs_ASY1_H2AZ_H3K27me3_peakOverlap_summits",
                   as.character(size), "bp/")
system(paste0("[ -d ", regionDir, " ] || mkdir ", regionDir))
system(paste0("[ -d ", motifDir, " ] || mkdir ", motifDir))

# peaks
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/",
                           markName, "/snakemake_ChIPseq/mapped/",
                           featureName, "peakProfiles_overlap_H2AZ_H3K27me3_peaks/heatmaps/gff/bodies_by_log2_",
                           libName, "_control/",
                           featureLibName, "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel, "_noMinWidth_in_",
                           genomeName, "genome_", region,
                           "_overlapping_H2AZ_Rep1_ChIP_and_H3K27me3_ChIP_SRR6350666_peaks",
                           "_ordered_by_log2_", libName, "_control_in_bodies.gff"))
colnames(peaks) <- c("chr", "source", "feature",
                     "start", "end", "sigval",
                     "strand", "qval", "summit0based")
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start,
                                    end = peaks$end),
                   strand = "*",
                   sigval = peaks$sigval,
                   qval = peaks$qval,
                   summit0based = peaks$summit0based)

# Extract peak summits +/- (size/2) bp for use in weeder2 motif analysis
peaksGR_summits <- GRanges(seqnames = seqnames(peaksGR),
                           ranges = IRanges(start = start(peaksGR)+peaksGR$summit0based,
                                            end = start(peaksGR)+peaksGR$summit0based),
                           strand = "*",
                           sigval = peaksGR$sigval,
                           qval = peaksGR$qval,
                           summit0based = peaksGR$summit0based)
peaksGR_summitsFlank <- locMidpointFlank(x = peaksGR_summits,
                                         leftFlank = (size/2),
                                         rightFlank = (size/2))
peaks_summitsFlank_bed <- data.frame(chr = seqnames(peaksGR_summitsFlank),
                                     start = start(peaksGR_summitsFlank)-1,
                                     end = end(peaksGR_summitsFlank),
                                     name = paste0(as.character(seqnames(peaksGR)),
                                                   ":",
                                                   as.character(start(peaksGR)),
                                                   "-",
                                                   as.character(end(peaksGR))))
write.table(peaks_summitsFlank_bed,
            file = paste0(motifDir, featureLibName, "_peaks_",
                          libName, "sorted_in_",
                          genomeName, "genome_", region,
                          "_summits", as.character(size), "bp.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")

# ranLoc
ranLoc <- read.table(paste0("/home/ajt200/analysis/wheat/",
                            markName, "/snakemake_ChIPseq/mapped/",
                            featureName, "peakProfiles_overlap_H2AZ_H3K27me3_peaks/heatmaps/gff/bodies_by_log2_",
                            libName, "_control/",
                            featureLibName, "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel, "_noMinWidth_in_",
                            genomeName, "genome_", region,
                            "_overlapping_H2AZ_Rep1_ChIP_and_H3K27me3_ChIP_SRR6350666_peaks_randomLoci.gff"))
colnames(ranLoc) <- c("chr", "source", "feature",
                      "start", "end", "score",
                      "strand", "phase", "attribute")
ranLocGR <- GRanges(seqnames = ranLoc$chr,
                    ranges = IRanges(start = ranLoc$start,
                                     end = ranLoc$end),
                    strand = "*")

# Extract peak summits +/- (size/2) bp for use in weeder2 motif analysis
ranLocGR_summits <- GRanges(seqnames = seqnames(ranLocGR),
                            ranges = IRanges(start = start(ranLocGR)+round((end(ranLocGR)-start(ranLocGR))/2),
                                             end = start(ranLocGR)+round((end(ranLocGR)-start(ranLocGR))/2)),
                            strand = "*")
ranLocGR_summitsFlank <- locMidpointFlank(x = ranLocGR_summits,
                                          leftFlank = (size/2),
                                          rightFlank = (size/2))
ranLoc_summitsFlank_bed <- data.frame(chr = seqnames(ranLocGR_summitsFlank),
                                      start = start(ranLocGR_summitsFlank)-1,
                                      end = end(ranLocGR_summitsFlank),
                                      name = paste0(as.character(seqnames(ranLocGR)),
                                                    ":",
                                                    as.character(start(ranLocGR)),
                                                    "-",
                                                    as.character(end(ranLocGR))))
write.table(ranLoc_summitsFlank_bed,
            file = paste0(motifDir, featureLibName, "_peaks_",
                          libName, "sorted_in_",
                          genomeName, "genome_", region,
                          "_summits", as.character(size), "bp_randomLoci.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")
