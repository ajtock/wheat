#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./varietal_SNP_frequency_chrProfiles.R 1Mb 1000000

#winName <- "1Mb"
#winSize <- 1000000

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])

library(rtracklayer)
library(parallel)
library(GenomicRanges)
library(doParallel)

SNPs <- read.table("all_filtered_SNPs_allaccessions_allploidy_SNPeff.vcf",
                   header = F, skip = 6,
                   colClasses = c(rep(NA, 2), rep("NULL", 67)))
colnames(SNPs) <- c("chr", "pos")

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
#chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
#chrLens <- chrLens[-length(chrLens)]

# Define windows as GRanges object
windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  seqWindows <- seq(1, chrLens[i], by = winSize)
  windowsIR <- IRanges(start = seqWindows,
                       width = winSize)
  windowsIR <- windowsIR[-length(windowsIR)]
  windowsIR <- append(windowsIR,
                      IRanges(start = seqWindows[length(seqWindows)],
                              end = chrLens[i]))
  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = windowsIR,
                          strand = "*")
  print(chrWindowsGR)
  windowsGR <- append(windowsGR, chrWindowsGR)
}

featureProfile <- NULL
for(i in 1:length(chrs)) {
  # Count SNPs within windows
  chrSNPs <- SNPs[SNPs$chr == chrs[i],]
  chrSNPsGR <- GRanges(seqnames = chrs[i],
                       ranges = IRanges(start = chrSNPs$pos,
                                        end = chrSNPs$pos),
                       strand = "*")
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  winSNPs <- countOverlaps(chrWindowsGR,
                           chrSNPsGR,
                           ignore.strand = T)
  chrProfile <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(chrWindowsGR)),
                           winSNPs = as.integer(winSNPs))
  featureProfile <- rbind(featureProfile, chrProfile)
}
write.table(featureProfile,
            file = paste0("varietal_SNP_frequency_per_", winName,
                          ".txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)
