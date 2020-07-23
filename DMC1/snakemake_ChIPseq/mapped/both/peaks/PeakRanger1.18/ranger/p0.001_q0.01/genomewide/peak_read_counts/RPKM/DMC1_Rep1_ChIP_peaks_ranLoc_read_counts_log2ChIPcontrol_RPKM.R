#!/applications/R/R-3.5.0/bin/Rscript

# Get read counts for each peak and ranLoc (with equivalent width distribution to peaks)
# and plot summaries

# Usage:
# ./DMC1_Rep1_ChIP_peaks_ranLoc_read_counts_log2ChIPcontrol_RPKM.R DMC1_Rep1_ChIP DMC1 'chr3B'

#libName <- "DMC1_Rep1_ChIP"
#dirName <- "DMC1"
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
                    R2C[x], ".gff"),
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

# Calcualte TPM (transcripts per kilobase per million) for each R1R3peak and R1R3ranLoc
R1R3peak_ChIP_RPK <- R1R3peak_ChIP_reads/(width(R1R3peaksGR)/1e+03)
R1R3peak_ChIP_RPKPM_scaling_factor <- sum(R1R3peak_ChIP_RPK)/1e+06
R1R3peak_ChIP_TPM <- R1R3peak_ChIP_RPK/R1R3peak_ChIP_RPKPM_scaling_factor

R1R3ranLoc_ChIP_RPK <- R1R3ranLoc_ChIP_reads/(width(R1R3ranLocGR)/1e+03)
R1R3ranLoc_ChIP_RPKPM_scaling_factor <- sum(R1R3ranLoc_ChIP_RPK)/1e+06
R1R3ranLoc_ChIP_TPM <- R1R3ranLoc_ChIP_RPK/R1R3ranLoc_ChIP_RPKPM_scaling_factor

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

# Calcualte TPM (transcripts per kilobase per million) for each R2Cpeak and R2CranLoc
R2Cpeak_ChIP_RPK <- R2Cpeak_ChIP_reads/(width(R2CpeaksGR)/1e+03)
R2Cpeak_ChIP_RPKPM_scaling_factor <- sum(R2Cpeak_ChIP_RPK)/1e+06
R2Cpeak_ChIP_TPM <- R2Cpeak_ChIP_RPK/R2Cpeak_ChIP_RPKPM_scaling_factor

R2CranLoc_ChIP_RPK <- R2CranLoc_ChIP_reads/(width(R2CranLocGR)/1e+03)
R2CranLoc_ChIP_RPKPM_scaling_factor <- sum(R2CranLoc_ChIP_RPK)/1e+06
R2CranLoc_ChIP_TPM <- R2CranLoc_ChIP_RPK/R2CranLoc_ChIP_RPKPM_scaling_factor


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

# Calcualte TPM (transcripts per kilobase per million) for each R1R3peak and R1R3ranLoc
R1R3peak_input_RPK <- R1R3peak_input_reads/(width(R1R3peaksGR)/1e+03)
R1R3peak_input_RPKPM_scaling_factor <- sum(R1R3peak_input_RPK)/1e+06
R1R3peak_input_TPM <- R1R3peak_input_RPK/R1R3peak_input_RPKPM_scaling_factor

R1R3ranLoc_input_RPK <- R1R3ranLoc_input_reads/(width(R1R3ranLocGR)/1e+03)
R1R3ranLoc_input_RPKPM_scaling_factor <- sum(R1R3ranLoc_input_RPK)/1e+06
R1R3ranLoc_input_TPM <- R1R3ranLoc_input_RPK/R1R3ranLoc_input_RPKPM_scaling_factor

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

# Calcualte TPM (transcripts per kilobase per million) for each R2Cpeak and R2CranLoc
R2Cpeak_input_RPK <- R2Cpeak_input_reads/(width(R2CpeaksGR)/1e+03)
R2Cpeak_input_RPKPM_scaling_factor <- sum(R2Cpeak_input_RPK)/1e+06
R2Cpeak_input_TPM <- R2Cpeak_input_RPK/R2Cpeak_input_RPKPM_scaling_factor

R2CranLoc_input_RPK <- R2CranLoc_input_reads/(width(R2CranLocGR)/1e+03)
R2CranLoc_input_RPKPM_scaling_factor <- sum(R2CranLoc_input_RPK)/1e+06
R2CranLoc_input_TPM <- R2CranLoc_input_RPK/R2CranLoc_input_RPKPM_scaling_factor



# Define function for making colours transparent (for peak width histograms)
makeTransparent <- function(thisColour, alpha = 210)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}



# Plot peak width histogram, and RPKM+1 or RPM+1 vs loci width and cumulative fraction of loci (peaks and random)
pdf(file = paste0(libName, "_",
                  paste0(featureName, collapse = "_"),
                  "_peak_width_hist_and_RPKMplus1_or_RPMplus1_vs_ecdf_and_locus_width_",
                  chrs, ".pdf"),
                  height = 8, width = 12)
par(mfrow = c(2, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
# RPKM+1
hist(width(R1R3peaksGR), breaks = 1000, col = makeTransparent("red"), border = NA, lwd = 2,
     xlab = "Peak width (bp)", ylab = "Peaks", main = "", cex.lab = 2, cex.axis = 2,
     xlim = c(0, 1500))
abline(v = mean(width(peaksGR)), col = "red", lty = 2, lwd = 2)
legend("topright",
       legend = paste0("Mean = ", mean(width(peaksGR))),
       col = "white",
       text.col = "red",
       ncol = 1, cex = 1, lwd = 2, bty = "n")

plot(ecdf(ranLoc_ChIP_RPKMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "black")
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(peak_ChIP_RPKMplus1), log = "x", xlim = c(1, 4), xaxt = "n", yaxt = "n", pch = ".",   
     xlab = "Coverage (RPKM+1)", ylab = "Cumulative fraction of loci", main = "", cex.lab = 2, col = "red")
legend("topright",
       legend = c("Peaks", "Random"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = ranLoc_ChIP_RPKMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 4), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = "black")
axis(side = 1, at = c(1:4), labels = c(1:4), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = peak_ChIP_RPKMplus1, y = width(peaksGR), pch = ".", log = "xy",
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

plot(ecdf(ranLoc_ChIP_RPMplus1), log = "x", xlim = c(1, 10), xaxt = "n", yaxt = "n", pch = ".",
     xlab = "", ylab = "", main = "", col = "black")
axis(side = 1, at = c(1:10), labels = c(1:10), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(peak_ChIP_RPMplus1), log = "x", xlim = c(1, 10), xaxt = "n", yaxt = "n", pch = ".",   
     xlab = "Coverage (RPM+1)", ylab = "Cumulative fraction of loci", main = "", cex.lab = 2, col = "red")
legend("topright",
       legend = c("Peaks", "Random"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 2, bty = "n")
box(lwd = 2)

plot(x = ranLoc_ChIP_RPMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 10), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = "black")
axis(side = 1, at = c(1:10), labels = c(1:10), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = peak_ChIP_RPMplus1, y = width(peaksGR), pch = ".", log = "xy",
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
