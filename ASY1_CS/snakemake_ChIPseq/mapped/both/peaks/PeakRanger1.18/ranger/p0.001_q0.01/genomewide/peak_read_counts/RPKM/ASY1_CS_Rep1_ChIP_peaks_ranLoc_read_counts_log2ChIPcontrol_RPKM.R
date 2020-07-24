#!/applications/R/R-3.5.0/bin/Rscript

# Get read counts for each peak and ranLoc (with equivalent width distribution to peaks)
# and plot summaries

# Usage:
# ./ASY1_CS_Rep1_ChIP_peaks_ranLoc_read_counts_log2ChIPcontrol_RPKM.R ASY1_CS_Rep1_ChIP ASY1_CS 'chr3B'

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#chrs <- unlist(strsplit("chr3B",
#                        split = ","))

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
chrs <- unlist(strsplit(args[3],
                        split = ","))

library(GenomicAlignments)
library(ShortRead)
library(rtracklayer)
library(regioneR)

inDir <- paste0("/home/ajt200/analysis/wheat/", dirName,
                "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/")
R1R3 <- unlist(strsplit(c("Agenome_distal,Bgenome_distal,Dgenome_distal"),
                        split = ","))
R2C <- unlist(strsplit(c("Agenome_interstitial,Bgenome_interstitial,Dgenome_interstitial,Agenome_proximal,Bgenome_proximal,Dgenome_proximal"),
                       split = ","))
# Load R1R3peaks
R1R3peaks <- lapply(seq_along(R1R3), function(x) {
  read.table(paste0(inDir,
                    libName,
                    "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                    R1R3[x], ".gff"),
             header = F, stringsAsFactors = F)
})
# Concatenate if R1R3peaks is a list of R1R3peak sets
if(length(R1R3peaks) > 1) {
  R1R3peaks <- do.call(rbind, R1R3peaks)
} else {
  R1R3peaks <- R1R3peaks[[1]]
}

# Convert R1R3peaks into GRanges
R1R3peaksGR <- GRanges(seqnames = R1R3peaks[,1],
                       ranges = IRanges(start = R1R3peaks[,4],
                                        end = R1R3peaks[,5]),
                       strand = "*")

# Load R2Cpeaks
R2Cpeaks <- lapply(seq_along(R2C), function(x) {
  read.table(paste0(inDir,
                    libName,
                    "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                    R2C[x], ".gff"),
             header = F, stringsAsFactors = F)
})
# Concatenate if R2Cpeaks is a list of R2Cpeak sets
if(length(R2Cpeaks) > 1) {
  R2Cpeaks <- do.call(rbind, R2Cpeaks)
} else {
  R2Cpeaks <- R2Cpeaks[[1]]
}

# Convert R2Cpeaks into GRanges
R2CpeaksGR <- GRanges(seqnames = R2Cpeaks[,1],
                      ranges = IRanges(start = R2Cpeaks[,4],
                                       end = R2Cpeaks[,5]),
                      strand = "*")

# Subset peaks to specified chromosome(s)
if(chrs != "allchrs") {
  R1R3peaksGR <- R1R3peaksGR[seqnames(R1R3peaksGR) %in% chrs]
  R2CpeaksGR <- R2CpeaksGR[seqnames(R2CpeaksGR) %in% chrs]
}

print(paste0(libName, " R1R3peaks width mean ="))
print(mean(width(R1R3peaksGR)))
print(paste0(libName, " R1R3peaks width range ="))
print(range(width(R1R3peaksGR)))

print(paste0(libName, " R2Cpeaks width mean ="))
print(mean(width(R2CpeaksGR)))
print(paste0(libName, " R2Cpeaks width range ="))
print(range(width(R2CpeaksGR)))


# Load R1R3ranLoc
R1R3ranLoc <- lapply(seq_along(R1R3), function(x) {
  read.table(paste0(inDir,
                    libName,
                    "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                    R1R3[x], "_randomLoci.gff"),
             header = F, stringsAsFactors = F)
})
# Concatenate if R1R3ranLoc is a list of R1R3peak sets
if(length(R1R3ranLoc) > 1) {
  R1R3ranLoc <- do.call(rbind, R1R3ranLoc)
} else {
  R1R3ranLoc <- R1R3ranLoc[[1]]
}

# Convert R1R3ranLoc into GRanges
R1R3ranLocGR <- GRanges(seqnames = R1R3ranLoc[,1],
                        ranges = IRanges(start = R1R3ranLoc[,4],
                                         end = R1R3ranLoc[,5]),
                        strand = "*")

# Load R2CranLoc
R2CranLoc <- lapply(seq_along(R2C), function(x) {
  read.table(paste0(inDir,
                    libName,
                    "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                    R2C[x], "_randomLoci.gff"),
             header = F, stringsAsFactors = F)
})
# Concatenate if R2CranLoc is a list of R2Cpeak sets
if(length(R2CranLoc) > 1) {
  R2CranLoc <- do.call(rbind, R2CranLoc)
} else {
  R2CranLoc <- R2CranLoc[[1]]
}

# Convert R2CranLoc into GRanges
R2CranLocGR <- GRanges(seqnames = R2CranLoc[,1],
                       ranges = IRanges(start = R2CranLoc[,4],
                                        end = R2CranLoc[,5]),
                       strand = "*")

# Subset ranLoc to specified chromosome(s)
if(chrs != "allchrs") {
  R1R3ranLocGR <- R1R3ranLocGR[seqnames(R1R3ranLocGR) %in% chrs]
  R2CranLocGR <- R2CranLocGR[seqnames(R2CranLocGR) %in% chrs]
}

print(paste0("Random loci for ", libName, " R1R3peaks width mean ="))
print(mean(width(R1R3ranLocGR)))
print(paste0("Random loci for ", libName, " R1R3peaks width range ="))
print(range(width(R1R3ranLocGR)))

print(paste0("Random loci for ", libName, " R2Cpeaks width mean ="))
print(mean(width(R2CranLocGR)))
print(paste0("Random loci for ", libName, " R2Cpeaks width range ="))
print(range(width(R2CranLocGR)))


# Load ChIP BAM file
#ChIP_reads <- readGAlignmentPairs(paste0("/home/ajt200/analysis/wheat/", dirName,
#                                    "/snakemake_ChIPseq/mapped/both/", libName,
#                                    "_MappedOn_wheat_v1.0_lowXM_both_sort.bam"))
#ChIP_reads <- ChIP_reads[seqnames(ChIP_reads) != "chrUn"]
#save(ChIP_reads,
#     file = paste0(libName, "_readGAlignmentPairs.RData"))
load(paste0(libName, "_readGAlignmentPairs.RData"))
ChIP_reads <- reads
rm(reads); gc()
ChIP_lib_size <- length(ChIP_reads)

# Calculate "per million" scaling factor
ChIP_RPM_scaling_factor <- ChIP_lib_size/1e6

if(chrs != "allchrs") {
  ChIP_reads <- ChIP_reads[seqnames(ChIP_reads) %in% chrs]
}

# Calculate RPM and RPKM for each R1R3peak and R1R3ranLoc
R1R3peak_ChIP_reads <- countOverlaps(query = R1R3peaksGR,
                                     subject = ChIP_reads)
R1R3peak_ChIP_RPM <- R1R3peak_ChIP_reads/ChIP_RPM_scaling_factor
R1R3peak_ChIP_RPKM <- R1R3peak_ChIP_RPM/(width(R1R3peaksGR)/1e+03)
R1R3peak_ChIP_RPMplus1 <- R1R3peak_ChIP_RPM+1
R1R3peak_ChIP_RPKMplus1 <- R1R3peak_ChIP_RPKM+1

R1R3ranLoc_ChIP_reads <- countOverlaps(query = R1R3ranLocGR,
                                       subject = ChIP_reads)
R1R3ranLoc_ChIP_RPM <- R1R3ranLoc_ChIP_reads/ChIP_RPM_scaling_factor
R1R3ranLoc_ChIP_RPKM <- R1R3ranLoc_ChIP_RPM/(width(R1R3ranLocGR)/1e+03)
R1R3ranLoc_ChIP_RPMplus1 <- R1R3ranLoc_ChIP_RPM+1
R1R3ranLoc_ChIP_RPKMplus1 <- R1R3ranLoc_ChIP_RPKM+1

# Calculate RPM and RPKM for each R2Cpeak and R2CranLoc
R2Cpeak_ChIP_reads <- countOverlaps(query = R2CpeaksGR,
                                    subject = ChIP_reads)
R2Cpeak_ChIP_RPM <- R2Cpeak_ChIP_reads/ChIP_RPM_scaling_factor
R2Cpeak_ChIP_RPKM <- R2Cpeak_ChIP_RPM/(width(R2CpeaksGR)/1e+03)
R2Cpeak_ChIP_RPMplus1 <- R2Cpeak_ChIP_RPM+1
R2Cpeak_ChIP_RPKMplus1 <- R2Cpeak_ChIP_RPKM+1

R2CranLoc_ChIP_reads <- countOverlaps(query = R2CranLocGR,
                                      subject = ChIP_reads)
R2CranLoc_ChIP_RPM <- R2CranLoc_ChIP_reads/ChIP_RPM_scaling_factor
R2CranLoc_ChIP_RPKM <- R2CranLoc_ChIP_RPM/(width(R2CranLocGR)/1e+03)
R2CranLoc_ChIP_RPMplus1 <- R2CranLoc_ChIP_RPM+1
R2CranLoc_ChIP_RPKMplus1 <- R2CranLoc_ChIP_RPKM+1


# Load input BAM file
#input_reads <- readGAlignmentPairs(paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/input/snakemake_ChIPseq/mapped/both/",
#                                          "H3_input_SRR6350669_MappedOn_wheat_v1.0_lowXM_both_sort.bam")) 
#input_reads <- input_reads[seqnames(input_reads) != "chrUn"]
#save(input_reads,
#     file = paste0("H3_input_SRR6350669_readGAlignmentPairs.RData"))
load(paste0("H3_input_SRR6350669_readGAlignmentPairs.RData"))
input_lib_size <- length(input_reads)

# Calculate "per million" scaling factor
input_RPM_scaling_factor <- input_lib_size/1e6

if(chrs != "allchrs") {
  input_reads <- input_reads[seqnames(input_reads) %in% chrs]
}

# Calculate RPM and RPKM for each R1R3peak and R1R3ranLoc
R1R3peak_input_reads <- countOverlaps(query = R1R3peaksGR,
                                      subject = input_reads)
R1R3peak_input_RPM <- R1R3peak_input_reads/input_RPM_scaling_factor
R1R3peak_input_RPKM <- R1R3peak_input_RPM/(width(R1R3peaksGR)/1e+03)
R1R3peak_input_RPMplus1 <- R1R3peak_input_RPM+1
R1R3peak_input_RPKMplus1 <- R1R3peak_input_RPKM+1

R1R3ranLoc_input_reads <- countOverlaps(query = R1R3ranLocGR,
                                        subject = input_reads)
R1R3ranLoc_input_RPM <- R1R3ranLoc_input_reads/input_RPM_scaling_factor
R1R3ranLoc_input_RPKM <- R1R3ranLoc_input_RPM/(width(R1R3ranLocGR)/1e+03)
R1R3ranLoc_input_RPMplus1 <- R1R3ranLoc_input_RPM+1
R1R3ranLoc_input_RPKMplus1 <- R1R3ranLoc_input_RPKM+1

# Calculate RPM and RPKM for each R2Cpeak and R2CranLoc
R2Cpeak_input_reads <- countOverlaps(query = R2CpeaksGR,
                                     subject = input_reads)
R2Cpeak_input_RPM <- R2Cpeak_input_reads/input_RPM_scaling_factor
R2Cpeak_input_RPKM <- R2Cpeak_input_RPM/(width(R2CpeaksGR)/1e+03)
R2Cpeak_input_RPMplus1 <- R2Cpeak_input_RPM+1
R2Cpeak_input_RPKMplus1 <- R2Cpeak_input_RPKM+1

R2CranLoc_input_reads <- countOverlaps(query = R2CranLocGR,
                                       subject = input_reads)
R2CranLoc_input_RPM <- R2CranLoc_input_reads/input_RPM_scaling_factor
R2CranLoc_input_RPKM <- R2CranLoc_input_RPM/(width(R2CranLocGR)/1e+03)
R2CranLoc_input_RPMplus1 <- R2CranLoc_input_RPM+1
R2CranLoc_input_RPKMplus1 <- R2CranLoc_input_RPKM+1

# Calculate log2(ChIP/input) coverage for R1R3peaks and R1R3ranLoc
log2_R1R3peak_ChIP_input_RPKMplus1 <- log2(R1R3peak_ChIP_RPKMplus1/R1R3peak_input_RPKMplus1)
log2_R1R3peak_ChIP_input_RPMplus1 <- log2(R1R3peak_ChIP_RPMplus1/R1R3peak_input_RPMplus1)

log2_R1R3ranLoc_ChIP_input_RPKMplus1 <- log2(R1R3ranLoc_ChIP_RPKMplus1/R1R3ranLoc_input_RPKMplus1)
log2_R1R3ranLoc_ChIP_input_RPMplus1 <- log2(R1R3ranLoc_ChIP_RPMplus1/R1R3ranLoc_input_RPMplus1)

# Calculate log2(ChIP/input) coverage for R2Cpeaks and R2CranLoc
log2_R2Cpeak_ChIP_input_RPKMplus1 <- log2(R2Cpeak_ChIP_RPKMplus1/R2Cpeak_input_RPKMplus1)
log2_R2Cpeak_ChIP_input_RPMplus1 <- log2(R2Cpeak_ChIP_RPMplus1/R2Cpeak_input_RPMplus1)

log2_R2CranLoc_ChIP_input_RPKMplus1 <- log2(R2CranLoc_ChIP_RPKMplus1/R2CranLoc_input_RPKMplus1)
log2_R2CranLoc_ChIP_input_RPMplus1 <- log2(R2CranLoc_ChIP_RPMplus1/R2CranLoc_input_RPMplus1)


# Define function for making colours transparent (for peak width histograms)
makeTransparent <- function(thisColour, alpha = 180)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}

# Calculate peak heat vs width correlation coefficients
# Standardise P-values to a sample size of 100 (q-values) as proposed by
# Good (1982) Standardized tail-area probabilities. Journal of Computation and Simulation 16: 65-66
# and summarised by Woolley (2003):
# https://stats.stackexchange.com/questions/22233/how-to-choose-significance-level-for-a-large-data-set
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.518.5341&rep=rep1&type=pdf
# Woolley (2003): "Clearly, the meaningfulness of the p-value diminishes as the sample size increases";
# Anne Z. (2012, Pearson eCollege, Denver): "In the real world, there are unlikely to be semi-partial correlations
# that are exactly zero, which is the null hypothesis in testing significance of a regression coefficient."
# Formally, the standardised p-value is defined as:
# q = min(0.5, p * sqrt( (n/100) ))
# Woolley (2003): "The value of 0.5 is somewhat arbitrary, though its purpose is to avoid q-values of greater than 1."
R1R3peak_L2FC_RPKMplus1_r <- round(cor.test(x = log2_R1R3peak_ChIP_input_RPKMplus1, y = width(R1R3peaksGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3peak_L2FC_RPKMplus1_p <- round(min(0.5, cor.test(x = log2_R1R3peak_ChIP_input_RPKMplus1, y = width(R1R3peaksGR), method = "spearman")$p.value * sqrt( (length(R1R3peaksGR)/100) )), digits = 4)
R2Cpeak_L2FC_RPKMplus1_r <- round(cor.test(x = log2_R2Cpeak_ChIP_input_RPKMplus1, y = width(R2CpeaksGR), method = "spearman")$estimate[[1]], digits = 2)
R2Cpeak_L2FC_RPKMplus1_p <- round(min(0.5, cor.test(x = log2_R2Cpeak_ChIP_input_RPKMplus1, y = width(R2CpeaksGR), method = "spearman")$p.value * sqrt( (length(R2CpeaksGR)/100) )), digits = 4)
R1R3ranLoc_L2FC_RPKMplus1_r <- round(cor.test(x = log2_R1R3ranLoc_ChIP_input_RPKMplus1, y = width(R1R3ranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3ranLoc_L2FC_RPKMplus1_p <- round(min(0.5, cor.test(x = log2_R1R3ranLoc_ChIP_input_RPKMplus1, y = width(R1R3ranLocGR), method = "spearman")$p.value * sqrt( (length(R1R3ranLocGR)/100) )), digits = 4)
R2CranLoc_L2FC_RPKMplus1_r <- round(cor.test(x = log2_R2CranLoc_ChIP_input_RPKMplus1, y = width(R2CranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R2CranLoc_L2FC_RPKMplus1_p <- round(min(0.5, cor.test(x = log2_R2CranLoc_ChIP_input_RPKMplus1, y = width(R2CranLocGR), method = "spearman")$p.value * sqrt( (length(R2CranLocGR)/100) )), digits = 4)

R1R3peak_L2FC_RPMplus1_r <- round(cor.test(x = log2_R1R3peak_ChIP_input_RPMplus1, y = width(R1R3peaksGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3peak_L2FC_RPMplus1_p <- round(min(0.5, cor.test(x = log2_R1R3peak_ChIP_input_RPMplus1, y = width(R1R3peaksGR), method = "spearman")$p.value * sqrt( (length(R1R3peaksGR)/100) )), digits = 4)
R2Cpeak_L2FC_RPMplus1_r <- round(cor.test(x = log2_R2Cpeak_ChIP_input_RPMplus1, y = width(R2CpeaksGR), method = "spearman")$estimate[[1]], digits = 2)
R2Cpeak_L2FC_RPMplus1_p <- round(min(0.5, cor.test(x = log2_R2Cpeak_ChIP_input_RPMplus1, y = width(R2CpeaksGR), method = "spearman")$p.value * sqrt( (length(R2CpeaksGR)/100) )), digits = 4)
R1R3ranLoc_L2FC_RPMplus1_r <- round(cor.test(x = log2_R1R3ranLoc_ChIP_input_RPMplus1, y = width(R1R3ranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3ranLoc_L2FC_RPMplus1_p <- round(min(0.5, cor.test(x = log2_R1R3ranLoc_ChIP_input_RPMplus1, y = width(R1R3ranLocGR), method = "spearman")$p.value * sqrt( (length(R1R3ranLocGR)/100) )), digits = 4)
R2CranLoc_L2FC_RPMplus1_r <- round(cor.test(x = log2_R2CranLoc_ChIP_input_RPMplus1, y = width(R2CranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R2CranLoc_L2FC_RPMplus1_p <- round(min(0.5, cor.test(x = log2_R2CranLoc_ChIP_input_RPMplus1, y = width(R2CranLocGR), method = "spearman")$p.value * sqrt( (length(R2CranLocGR)/100) )), digits = 4)

R1R3peak_ChIP_RPKMplus1_r <- round(cor.test(x = R1R3peak_ChIP_RPKMplus1, y = width(R1R3peaksGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3peak_ChIP_RPKMplus1_p <- round(min(0.5, cor.test(x = R1R3peak_ChIP_RPKMplus1, y = width(R1R3peaksGR), method = "spearman")$p.value * sqrt( (length(R1R3peaksGR)/100) )), digits = 4)
R2Cpeak_ChIP_RPKMplus1_r <- round(cor.test(x = R2Cpeak_ChIP_RPKMplus1, y = width(R2CpeaksGR), method = "spearman")$estimate[[1]], digits = 2)
R2Cpeak_ChIP_RPKMplus1_p <- round(min(0.5, cor.test(x = R2Cpeak_ChIP_RPKMplus1, y = width(R2CpeaksGR), method = "spearman")$p.value * sqrt( (length(R2CpeaksGR)/100) )), digits = 4)
R1R3ranLoc_ChIP_RPKMplus1_r <- round(cor.test(x = R1R3ranLoc_ChIP_RPKMplus1, y = width(R1R3ranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3ranLoc_ChIP_RPKMplus1_p <- round(min(0.5, cor.test(x = R1R3ranLoc_ChIP_RPKMplus1, y = width(R1R3ranLocGR), method = "spearman")$p.value * sqrt( (length(R1R3ranLocGR)/100) )), digits = 4)
R2CranLoc_ChIP_RPKMplus1_r <- round(cor.test(x = R2CranLoc_ChIP_RPKMplus1, y = width(R2CranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R2CranLoc_ChIP_RPKMplus1_p <- round(min(0.5, cor.test(x = R2CranLoc_ChIP_RPKMplus1, y = width(R2CranLocGR), method = "spearman")$p.value * sqrt( (length(R2CranLocGR)/100) )), digits = 4)

R1R3peak_ChIP_RPMplus1_r <- round(cor.test(x = R1R3peak_ChIP_RPMplus1, y = width(R1R3peaksGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3peak_ChIP_RPMplus1_p <- round(min(0.5, cor.test(x = R1R3peak_ChIP_RPMplus1, y = width(R1R3peaksGR), method = "spearman")$p.value * sqrt( (length(R1R3peaksGR)/100) )), digits = 4)
R2Cpeak_ChIP_RPMplus1_r <- round(cor.test(x = R2Cpeak_ChIP_RPMplus1, y = width(R2CpeaksGR), method = "spearman")$estimate[[1]], digits = 2)
R2Cpeak_ChIP_RPMplus1_p <- round(min(0.5, cor.test(x = R2Cpeak_ChIP_RPMplus1, y = width(R2CpeaksGR), method = "spearman")$p.value * sqrt( (length(R2CpeaksGR)/100) )), digits = 4)
R1R3ranLoc_ChIP_RPMplus1_r <- round(cor.test(x = R1R3ranLoc_ChIP_RPMplus1, y = width(R1R3ranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R1R3ranLoc_ChIP_RPMplus1_p <- round(min(0.5, cor.test(x = R1R3ranLoc_ChIP_RPMplus1, y = width(R1R3ranLocGR), method = "spearman")$p.value * sqrt( (length(R1R3ranLocGR)/100) )), digits = 4)
R2CranLoc_ChIP_RPMplus1_r <- round(cor.test(x = R2CranLoc_ChIP_RPMplus1, y = width(R2CranLocGR), method = "spearman")$estimate[[1]], digits = 2)
R2CranLoc_ChIP_RPMplus1_p <- round(min(0.5, cor.test(x = R2CranLoc_ChIP_RPMplus1, y = width(R2CranLocGR), method = "spearman")$p.value * sqrt( (length(R2CranLocGR)/100) )), digits = 4)


# Plot peak width histogram, and RPKM+1 or RPM+1 vs loci width and cumulative fraction of loci (peaks and random)
pdf(file = paste0(libName,
                  "_peak_width_hist_and_log2ChIPinput_RPKMplus1_or_RPMplus1_vs_ecdf_and_locus_width_",
                  chrs, ".pdf"),
                  height = 16, width = 12)
par(mfrow = c(4, 3), mar =  c(6, 6, 2, 2), mgp = c(3, 1, 0))
# log2(ChIP/input) RPKM
minBreak <- min(c(width(R1R3peaksGR), width(R2CpeaksGR))) - 0.001
maxBreak <- max(c(width(R1R3peaksGR), width(R2CpeaksGR)))
vecBreak <- pretty(minBreak:maxBreak, n = 1000)
histR1R3 <- hist(width(R1R3peaksGR), breaks = vecBreak, plot = F)
histR2C <- hist(width(R2CpeaksGR), breaks = vecBreak, plot = F)
plot(histR2C, col = makeTransparent("red4"), border = NA, lwd = 2,
     xlab = "Peak width (bp)", ylab = "Peaks", main = "", xlim = c(0, 1500),
     cex.lab = 2, cex.axis = 2)
plot(add = T, histR1R3, col = makeTransparent("red"), border = NA, lwd = 2, xlim = c(0, 1500))
abline(v = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("right",
       legend = c(paste0("R1 & R3 peaks mean = ", round(mean(width(R1R3peaksGR))), " bp"),
                  paste0("R2a-R2b peaks mean = ", round(mean(width(R2CpeaksGR))), " bp")),
       col = "white",
       text.col = c("red", "red4"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")

plot(ecdf(log2_R2CranLoc_ChIP_input_RPKMplus1), xlim = c(-1, 1.5), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "navy")
axis(side = 1, lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(log2_R1R3ranLoc_ChIP_input_RPKMplus1), xlim = c(-1, 1.5),xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "blue")
par(new = T)
plot(ecdf(log2_R2Cpeak_ChIP_input_RPKMplus1), xlim = c(-1, 1.5), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "red4")
par(new = T)
plot(ecdf(log2_R1R3peak_ChIP_input_RPKMplus1), xlim = c(-1, 1.5), xaxt = "n", yaxt = "n", pch = ".",
     xlab = expression("Log"[2]*"(ChIP/control) RPKM"),
     ylab = "Cumulative fraction of loci",
     main = "", cex.lab = 2, col = "red")
legend("right",
       legend = c("R1 & R3 peaks", "R2a-R2b peaks",
                  "R1 & R3 random loci", "R2a-R2b random loci"),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = log2_R2CranLoc_ChIP_input_RPKMplus1, y = width(R2CranLocGR), pch = ".", log = "y",
     xlim = c(-1, 1.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("navy"))
axis(side = 1, lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = log2_R1R3ranLoc_ChIP_input_RPKMplus1, y = width(R1R3ranLocGR), pch = ".", log = "y",
     xlim = c(-1, 1.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("blue"))
par(new = T)
plot(x = log2_R2Cpeak_ChIP_input_RPKMplus1, y = width(R2CpeaksGR), pch = ".", log = "y",
     xlim = c(-1, 1.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("red4"))
par(new = T)
plot(x = log2_R1R3peak_ChIP_input_RPKMplus1, y = width(R1R3peaksGR), pch = ".", log = "y",
     xlim = c(-1, 1.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = expression("Log"[2]*"(ChIP/control) RPKM"),
     ylab = "Locus width (bp)",
     cex.lab = 2, col = makeTransparent("red"))
abline(h = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("bottomright",
       legend = c(as.expression(bquote("R1 & R3 peaks" ~ italic("r"[s]) ~ "=" ~ .(R1R3peak_L2FC_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R1R3peak_L2FC_RPKMplus1_p))),
                  bquote("R2a-R2b peaks" ~ italic("r"[s]) ~ "=" ~ .(R2Cpeak_L2FC_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R2Cpeak_L2FC_RPKMplus1_p)),
                  bquote("R1 & R3 random loci" ~ italic("r"[s]) ~ "=" ~ .(R1R3ranLoc_L2FC_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R1R3ranLoc_L2FC_RPKMplus1_p)),
                  bquote("R2a-R2b random loci" ~ italic("r"[s]) ~ "=" ~ .(R2CranLoc_L2FC_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R2CranLoc_L2FC_RPKMplus1_p))),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

# log2(ChIP/input) RPM
minBreak <- min(c(width(R1R3peaksGR), width(R2CpeaksGR))) - 0.001
maxBreak <- max(c(width(R1R3peaksGR), width(R2CpeaksGR)))
vecBreak <- pretty(minBreak:maxBreak, n = 1000)
histR1R3 <- hist(width(R1R3peaksGR), breaks = vecBreak, plot = F)
histR2C <- hist(width(R2CpeaksGR), breaks = vecBreak, plot = F)
plot(histR2C, col = makeTransparent("red4"), border = NA, lwd = 2,
     xlab = "Peak width (bp)", ylab = "Peaks", main = "", xlim = c(0, 1500),
     cex.lab = 2, cex.axis = 2)
plot(add = T, histR1R3, col = makeTransparent("red"), border = NA, lwd = 2, xlim = c(0, 1500))
abline(v = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("right",
       legend = c(paste0("R1 & R3 peaks mean = ", round(mean(width(R1R3peaksGR))), " bp"),
                  paste0("R2a-R2b peaks mean = ", round(mean(width(R2CpeaksGR))), " bp")),
       col = "white",
       text.col = c("red", "red4"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")

plot(ecdf(log2_R2CranLoc_ChIP_input_RPMplus1), xlim = c(-0.5,0.5), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "navy")
axis(side = 1, lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(log2_R1R3ranLoc_ChIP_input_RPMplus1), xlim = c(-0.5,0.5),xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "blue")
par(new = T)
plot(ecdf(log2_R2Cpeak_ChIP_input_RPMplus1), xlim = c(-0.5,0.5), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "red4")
par(new = T)
plot(ecdf(log2_R1R3peak_ChIP_input_RPMplus1), xlim = c(-0.5,0.5), xaxt = "n", yaxt = "n", pch = ".",
     xlab = expression("Log"[2]*"(ChIP/control) RPM"),
     ylab = "Cumulative fraction of loci",
     main = "", cex.lab = 2, col = "red")
legend("right",
       legend = c("R1 & R3 peaks", "R2a-R2b peaks",
                  "R1 & R3 random loci", "R2a-R2b random loci"),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = log2_R2CranLoc_ChIP_input_RPMplus1, y = width(R2CranLocGR), pch = ".", log = "y",
     xlim = c(-0.5,0.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("navy"))
axis(side = 1, lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = log2_R1R3ranLoc_ChIP_input_RPMplus1, y = width(R1R3ranLocGR), pch = ".", log = "y",
     xlim = c(-0.5,0.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("blue"))
par(new = T)
plot(x = log2_R2Cpeak_ChIP_input_RPMplus1, y = width(R2CpeaksGR), pch = ".", log = "y",
     xlim = c(-0.5,0.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("red4"))
par(new = T)
plot(x = log2_R1R3peak_ChIP_input_RPMplus1, y = width(R1R3peaksGR), pch = ".", log = "y",
     xlim = c(-0.5,0.5), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = expression("Log"[2]*"(ChIP/control) RPM"),
     ylab = "Locus width (bp)",
     cex.lab = 2, col = makeTransparent("red"))
abline(h = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("bottomright",
       legend = c(as.expression(bquote("R1 & R3 peaks" ~ italic("r"[s]) ~ "=" ~ .(R1R3peak_L2FC_RPMplus1_r) * ";" ~ italic("P =") ~ .(R1R3peak_L2FC_RPMplus1_p))),
                  bquote("R2a-R2b peaks" ~ italic("r"[s]) ~ "=" ~ .(R2Cpeak_L2FC_RPMplus1_r) * ";" ~ italic("P =") ~ .(R2Cpeak_L2FC_RPMplus1_p)),
                  bquote("R1 & R3 random loci" ~ italic("r"[s]) ~ "=" ~ .(R1R3ranLoc_L2FC_RPMplus1_r) * ";" ~ italic("P =") ~ .(R1R3ranLoc_L2FC_RPMplus1_p)),
                  bquote("R2a-R2b random loci" ~ italic("r"[s]) ~ "=" ~ .(R2CranLoc_L2FC_RPMplus1_r) * ";" ~ italic("P =") ~ .(R2CranLoc_L2FC_RPMplus1_p))),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

# ChIP RPKM+1
minBreak <- min(c(width(R1R3peaksGR), width(R2CpeaksGR))) - 0.001
maxBreak <- max(c(width(R1R3peaksGR), width(R2CpeaksGR)))
vecBreak <- pretty(minBreak:maxBreak, n = 1000)
histR1R3 <- hist(width(R1R3peaksGR), breaks = vecBreak, plot = F)
histR2C <- hist(width(R2CpeaksGR), breaks = vecBreak, plot = F)
plot(histR2C, col = makeTransparent("red4"), border = NA, lwd = 2,
     xlab = "Peak width (bp)", ylab = "Peaks", main = "", xlim = c(0, 1500),
     cex.lab = 2, cex.axis = 2)
plot(add = T, histR1R3, col = makeTransparent("red"), border = NA, lwd = 2, xlim = c(0, 1500))
abline(v = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("right",
       legend = c(paste0("R1 & R3 peaks mean = ", round(mean(width(R1R3peaksGR))), " bp"),
                  paste0("R2a-R2b peaks mean = ", round(mean(width(R2CpeaksGR))), " bp")),
       col = "white",
       text.col = c("red", "red4"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")

plot(ecdf(R2CranLoc_ChIP_RPKMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "navy")
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(R1R3ranLoc_ChIP_RPKMplus1), log = "x", xlim = c(1, 4),xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "blue")
par(new = T)
plot(ecdf(R2Cpeak_ChIP_RPKMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "red4")
par(new = T)
plot(ecdf(R1R3peak_ChIP_RPKMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "Coverage RPKM+1",
     ylab = "Cumulative fraction of loci",
     main = "", cex.lab = 2, col = "red")
legend("right",
       legend = c("R1 & R3 peaks", "R2a-R2b peaks",
                  "R1 & R3 random loci", "R2a-R2b random loci"),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = R2CranLoc_ChIP_RPKMplus1, y = width(R2CranLocGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("navy"))
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = R1R3ranLoc_ChIP_RPKMplus1, y = width(R1R3ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("blue"))
par(new = T)
plot(x = R2Cpeak_ChIP_RPKMplus1, y = width(R2CpeaksGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("red4"))
par(new = T)
plot(x = R1R3peak_ChIP_RPKMplus1, y = width(R1R3peaksGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "Coverage RPKM+1",
     ylab = "Locus width (bp)",
     cex.lab = 2, col = makeTransparent("red"))
abline(h = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("bottomright",
       legend = c(as.expression(bquote("R1 & R3 peaks" ~ italic("r"[s]) ~ "=" ~ .(R1R3peak_ChIP_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R1R3peak_ChIP_RPKMplus1_p))),
                  bquote("R2a-R2b peaks" ~ italic("r"[s]) ~ "=" ~ .(R2Cpeak_ChIP_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R2Cpeak_ChIP_RPKMplus1_p)),
                  bquote("R1 & R3 random loci" ~ italic("r"[s]) ~ "=" ~ .(R1R3ranLoc_ChIP_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R1R3ranLoc_ChIP_RPKMplus1_p)),
                  bquote("R2a-R2b random loci" ~ italic("r"[s]) ~ "=" ~ .(R2CranLoc_ChIP_RPKMplus1_r) * ";" ~ italic("P =") ~ .(R2CranLoc_ChIP_RPKMplus1_p))),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

# ChIP RPM+1
minBreak <- min(c(width(R1R3peaksGR), width(R2CpeaksGR))) - 0.001
maxBreak <- max(c(width(R1R3peaksGR), width(R2CpeaksGR)))
vecBreak <- pretty(minBreak:maxBreak, n = 1000)
histR1R3 <- hist(width(R1R3peaksGR), breaks = vecBreak, plot = F)
histR2C <- hist(width(R2CpeaksGR), breaks = vecBreak, plot = F)
plot(histR2C, col = makeTransparent("red4"), border = NA, lwd = 2,
     xlab = "Peak width (bp)", ylab = "Peaks", main = "", xlim = c(0, 1500),
     cex.lab = 2, cex.axis = 2)
plot(add = T, histR1R3, col = makeTransparent("red"), border = NA, lwd = 2, xlim = c(0, 1500))
abline(v = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("right",
       legend = c(paste0("R1 & R3 peaks mean = ", round(mean(width(R1R3peaksGR))), " bp"),
                  paste0("R2a-R2b peaks mean = ", round(mean(width(R2CpeaksGR))), " bp")),
       col = "white",
       text.col = c("red", "red4"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")

plot(ecdf(R2CranLoc_ChIP_RPMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "navy")
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(R1R3ranLoc_ChIP_RPMplus1), log = "x", xlim = c(1, 4),xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "blue")
par(new = T)
plot(ecdf(R2Cpeak_ChIP_RPMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "red4")
par(new = T)
plot(ecdf(R1R3peak_ChIP_RPMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "Coverage RPM+1",
     ylab = "Cumulative fraction of loci",
     main = "", cex.lab = 2, col = "red")
legend("right",
       legend = c("R1 & R3 peaks", "R2a-R2b peaks",
                  "R1 & R3 random loci", "R2a-R2b random loci"),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = R2CranLoc_ChIP_RPMplus1, y = width(R2CranLocGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("navy"))
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = R1R3ranLoc_ChIP_RPMplus1, y = width(R1R3ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("blue"))
par(new = T)
plot(x = R2Cpeak_ChIP_RPMplus1, y = width(R2CpeaksGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = makeTransparent("red4"))
par(new = T)
plot(x = R1R3peak_ChIP_RPMplus1, y = width(R1R3peaksGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "Coverage RPM+1",
     ylab = "Locus width (bp)",
     cex.lab = 2, col = makeTransparent("red"))
abline(h = c(mean(width(R2CpeaksGR)), mean(width(R1R3peaksGR))), col = c("red4", "red"), lty = 2, lwd = 2)
legend("bottomright",
       legend = c(as.expression(bquote("R1 & R3 peaks" ~ italic("r"[s]) ~ "=" ~ .(R1R3peak_ChIP_RPMplus1_r) * ";" ~ italic("P =") ~ .(R1R3peak_ChIP_RPMplus1_p))),
                  bquote("R2a-R2b peaks" ~ italic("r"[s]) ~ "=" ~ .(R2Cpeak_ChIP_RPMplus1_r) * ";" ~ italic("P =") ~ .(R2Cpeak_ChIP_RPMplus1_p)),
                  bquote("R1 & R3 random loci" ~ italic("r"[s]) ~ "=" ~ .(R1R3ranLoc_ChIP_RPMplus1_r) * ";" ~ italic("P =") ~ .(R1R3ranLoc_ChIP_RPMplus1_p)),
                  bquote("R2a-R2b random loci" ~ italic("r"[s]) ~ "=" ~ .(R2CranLoc_ChIP_RPMplus1_r) * ";" ~ italic("P =") ~ .(R2CranLoc_ChIP_RPMplus1_p))),
       col = "white",
       text.col = c("red", "red4",
                    "blue", "navy"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

dev.off()
