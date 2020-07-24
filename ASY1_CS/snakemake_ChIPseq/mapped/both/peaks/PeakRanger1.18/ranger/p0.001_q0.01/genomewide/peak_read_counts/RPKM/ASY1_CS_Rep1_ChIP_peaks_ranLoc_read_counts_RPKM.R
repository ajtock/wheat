#!/applications/R/R-3.5.0/bin/Rscript

# Get read counts for each peak and ranLoc (with equivalent width distribution to peaks)
# and plot summaries

# Usage:
# ./ASY1_CS_Rep1_ChIP_peaks_ranLoc_read_counts_RPKM.R ASY1_CS_Rep1_ChIP ASY1_CS 'Agenome_distal,Bgenome_distal,Dgenome_distal,Agenome_interstitial,Bgenome_interstitial,Dgenome_interstitial,Agenome_proximal,Bgenome_proximal,Dgenome_proximal' 'chr3B'

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#featureName <- unlist(strsplit("Agenome_distal,Bgenome_distal,Dgenome_distal,Agenome_interstitial,Bgenome_interstitial,Dgenome_interstitial,Agenome_proximal,Bgenome_proximal,Dgenome_proximal",
#                               split = ","))
#chrs <- unlist(strsplit("chr3B",
#                        split = ","))

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
chrs <- unlist(strsplit(args[4],
                        split = ","))

library(GenomicAlignments)
library(ShortRead)
library(rtracklayer)
library(regioneR)


inDir <- paste0("/home/ajt200/analysis/wheat/", dirName,
                "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/")
# Load peaks
peaks <- lapply(seq_along(featureName), function(x) {
  read.table(paste0(inDir,
                    libName,
                    "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                    featureName[x], ".gff"),
             header = F, stringsAsFactors = F)
})
# Concatenate if peaks is a list of peak sets
if(length(peaks) > 1) {
  peaks <- do.call(rbind, peaks)
} else {
  peaks <- peaks[[1]]
}

# Convert peaks into GRanges
peaksGR <- GRanges(seqnames = peaks[,1],
                   ranges = IRanges(start = peaks[,4],
                                    end = peaks[,5]),
                   strand = "*")

# Subset peaks to specified chromosome(s)
if(chrs != "allchrs") {
  peaksGR <- peaksGR[seqnames(peaksGR) %in% chrs]
}

print(paste0(libName, " peaks width mean ="))
print(mean(width(peaksGR)))
print(paste0(libName, " peaks width range ="))
print(range(width(peaksGR)))

# Load ranLoc
ranLoc <- lapply(seq_along(featureName), function(x) {
  read.table(paste0(inDir,
                    libName,
                    "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                    featureName[x], "_randomLoci.gff"),
             header = F, stringsAsFactors = F)
})
# Concatenate if ranLoc is a list of peak sets
if(length(ranLoc) > 1) {
  ranLoc <- do.call(rbind, ranLoc)
} else {
  ranLoc <- ranLoc[[1]]
}

# Convert ranLoc into GRanges
ranLocGR <- GRanges(seqnames = ranLoc[,1],
                    ranges = IRanges(start = ranLoc[,4],
                                     end = ranLoc[,5]),
                    strand = "*")

# Subset ranLoc to specified chromosome(s)
if(chrs != "allchrs") {
  ranLocGR <- ranLocGR[seqnames(ranLocGR) %in% chrs]
}

print(paste0("Random loci for ", libName, " peaks width mean ="))
print(mean(width(ranLocGR)))
print(paste0("Random loci for ", libName, " peaks width range ="))
print(range(width(ranLocGR)))

# Load BAM file
#reads <- readGAlignmentPairs(paste0("/home/ajt200/analysis/wheat/", dirName,
#                                    "/snakemake_ChIPseq/mapped/both/", libName,
#                                    "_MappedOn_wheat_v1.0_lowXM_both_sort.bam"))
#reads <- reads[seqnames(reads) != "chrUn"]
#save(reads,
#     file = paste0(libName, "_readGAlignmentPairs.RData"))
load(paste0(libName, "_readGAlignmentPairs.RData"))
lib_size <- length(reads)

# Calculate "per million" scaling factor
RPM_scaling_factor <- lib_size/1e6

# Calculate RPM and RPKM for each peak and ranLoc
if(chrs != "allchrs") {
  reads <- reads[seqnames(reads) %in% chrs]
}
peak_reads <- countOverlaps(query = peaksGR,
                            subject = reads)
peak_RPM <- peak_reads/RPM_scaling_factor
peak_RPKM <- peak_RPM/(width(peaksGR)/1e+03)
peak_RPMplus1 <- peak_RPM+1
peak_RPKMplus1 <- peak_RPKM+1

ranLoc_reads <- countOverlaps(query = ranLocGR,
                              subject = reads)
ranLoc_RPM <- ranLoc_reads/RPM_scaling_factor
ranLoc_RPKM <- ranLoc_RPM/(width(ranLocGR)/1e+03)
ranLoc_RPMplus1 <- ranLoc_RPM+1
ranLoc_RPKMplus1 <- ranLoc_RPKM+1

# Calcualte TPM (transcripts per kilobase per million) for each peak and ranLoc
peak_RPK <- peak_reads/(width(peaksGR)/1e+03)
RPKPM_scaling_factor <- sum(peak_RPK)/1e+06
peak_TPM <- peak_RPK/RPKPM_scaling_factor

ranLoc_RPK <- ranLoc_reads/(width(ranLocGR)/1e+03)
RPKPM_scaling_factor <- sum(ranLoc_RPK)/1e+06
ranLoc_TPM <- ranLoc_RPK/RPKPM_scaling_factor


# Plot peak width histogram, and RPKM+1 or RPM+1 vs loci width and cumulative fraction of loci (peaks and random)
pdf(file = paste0(libName, "_",
                  paste0(featureName, collapse = "_"),
                  "_peak_width_hist_and_RPKMplus1_or_RPMplus1_vs_ecdf_and_locus_width_",
                  chrs, ".pdf"),
                  height = 8, width = 12)
par(mfrow = c(2, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
# RPKM+1
hist(width(peaksGR), breaks = 1000, col = "grey60", border = NA, lwd = 2,
     xlab = "Peak width (bp)", ylab = "Peaks", main = "", cex.lab = 2, cex.axis = 2,
     xlim = c(0, 1500))
abline(v = mean(width(peaksGR)), col = "red", lty = 2, lwd = 2)
legend("topright",
       legend = paste0("Mean = ", mean(width(peaksGR))),
       col = "white",
       text.col = "red",
       ncol = 1, cex = 1, lwd = 2, bty = "n")

plot(ecdf(ranLoc_RPKMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "black")
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(peak_RPKMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",   
     xlab = "Coverage (RPKM+1)", ylab = "Cumulative fraction of loci", main = "", cex.lab = 2, col = "red")
legend("topright",
       legend = c("Peaks", "Random"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = ranLoc_RPKMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = "black")
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = peak_RPKMplus1, y = width(peaksGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "Coverage (RPKM+1)", ylab = "Locus width (bp)", cex.lab = 2, col = "red")
abline(h = mean(width(peaksGR)), col = "grey60", lty = 2, lwd = 2)
legend("topright",
       legend = c("Peaks", "Random"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

# RPM+1
hist(width(peaksGR), breaks = 1000, col = "grey60", border = NA, lwd = 2,
     xlab = "Peak width (bp)", ylab = "Peaks", main = "", cex.lab = 2, cex.axis = 2,
     xlim = c(0, 1500))
abline(v = mean(width(peaksGR)), col = "red", lty = 2, lwd = 2)
legend("topright",
       legend = paste0("Mean = ", mean(width(peaksGR))),
       col = "white",
       text.col = "red",
       ncol = 1, cex = 1, lwd = 2, bty = "n")

plot(ecdf(ranLoc_RPMplus1), log = "x", xlim = c(1, 10), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "black")
axis(side = 1, at = c(1:10), labels = c(1:10), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(peak_RPMplus1), log = "x", xlim = c(1, 10), xaxt = "n", yaxt = "n", pch = ".",   
     xlab = "Coverage (RPM+1)", ylab = "Cumulative fraction of loci", main = "", cex.lab = 2, col = "red")
legend("topright",
       legend = c("Peaks", "Random"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = ranLoc_RPMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 10), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = "black")
axis(side = 1, at = c(1:10), labels = c(1:10), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = peak_RPMplus1, y = width(peaksGR), pch = ".", log = "xy",
     xlim = c(1, 10), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "Coverage (RPM+1)", ylab = "Locus width (bp)", cex.lab = 2, col = "red")
abline(h = mean(width(peaksGR)), col = "grey60", lty = 2, lwd = 2)
legend("topright",
       legend = c("Peaks", "Random"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

dev.off()
