#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./NLR_gene_frequency_chrProfiles.R 1Mb 1000000 15

#winName <- "1Mb"
#winSize <- 1000000
#N <- 15

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
N <- as.numeric(args[3])

library(rtracklayer)
library(parallel)
library(GenomicRanges)

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
outDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/genes/"
genesGFF <- readGFF(paste0(inDir, "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLR_genes_complete_representative_mRNA.gff3"))

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
  # Count genes within windows
  chrgenes <- genesGFF[genesGFF$seqid == chrs[i],]
  chrgenesGR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chrgenes$start,
                                         end = chrgenes$end),
                        strand = chrgenes$strand)
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  wingenes <- countOverlaps(chrWindowsGR,
                            chrgenesGR,
                            ignore.strand = T)
  chrProfile <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(chrWindowsGR)),
                           features = as.integer(wingenes))
  featureProfile <- rbind(featureProfile, chrProfile)
}
write.table(featureProfile,
            file = paste0(outDir,
                          "NLR_gene_frequency_per_", winName,
                          ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  featureProfile[featureProfile$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfiles <- mclapply(seq_along(chrProfiles), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfiles[[x]]$features,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfiles[[x]]$chr),
             window = as.integer(chrProfiles[[x]]$window),
             filt_features = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfiles))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_featureProfile <- do.call(rbind, filt_chrProfiles)

write.table(filt_featureProfile,
            file = paste0(outDir,
                          "NLR_gene_frequency_per_", winName,
                          "_smoothed.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
