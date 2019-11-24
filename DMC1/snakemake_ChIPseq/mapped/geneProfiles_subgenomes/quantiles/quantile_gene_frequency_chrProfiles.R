
#!/applications/R/R-3.5.0/bin/Rscript

# Calculate windowed feature frequencies for features within
# libName quantiles

# Usage:
# Rscript ./quantile_gene_frequency_chrProfiles.R ASY1_CS_Rep1_ChIP ASY1_CS 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' promoters 4 1Mb 1000000

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#region <- "promoters"
#quantiles <- 4
#winName <- "1Mb"
#winSize <- 1000000

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
featureName <- unlist(strsplit(args[3],
                               split = ","))
region <- args[4]
quantiles <- as.numeric(args[5])
winName <- args[6]
winSize <- as.numeric(args[7])

library(rtracklayer)
library(parallel)
library(GenomicRanges)

outDir <- paste0("quantiles_by_log2_", libName,
                 "_control_in_", region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

# Define plot titles
featureNamePlot <- paste0(sub("_\\w+", "", dirName), " ",
                          substr(featureName[1], start = 1, stop = 4),
                          " quantiles")
# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")

# Load table of features grouped into quantiles
# by decreasing log2(libName/control)
featuresDF <- read.table(paste0(outDir, "features_", quantiles, "quantiles",
                                "_by_log2_", libName, "_control_in_",
                                region, "_of_",
                                substring(featureName[1][1], first = 1, last = 5), "_in_",
                                paste0(substring(featureName, first = 10, last = 16),
                                       collapse = "_"), "_",
                                substring(featureName[1][1], first = 18), ".txt"),
                         header = T, sep = "\t")

quantilesDFlist <- mclapply(1:quantiles, function(k) {
  featuresDF[featuresDF$quantile == paste0("Quantile ", k),]
}, mc.cores = quantiles)

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

library(doParallel)
registerDoParallel(cores = quantiles)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(k = 1:quantiles) %dopar% {
  quantilekProfiles <- data.frame()
  for(i in 1:length(chrs)) {
    # Count quantile features within adjacent windows
    featuresChr <- quantilesDFlist[[k]][quantilesDFlist[[k]]$seqnames == chrs[i],]
    featuresChrGR <- GRanges(seqnames = chrs[i],
                             ranges = IRanges(start = featuresChr$start,
                                              end = featuresChr$end),
                             strand = featuresChr$strand)
    windowsChrGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
    winfeatures <- countOverlaps(query = windowsChrGR,
                                 subject = featuresChrGR,
                                 type = "any",
                                 ignore.strand = T)
    profileChr <- data.frame(chr = as.character(chrs[i]),
                             window = as.integer(start(windowsChrGR)),
                             winfeatures = as.integer(winfeatures))
    quantilekProfiles <- rbind(quantilekProfiles, profileChr)
  }
  write.table(quantilekProfiles,
              file = paste0(outDir,
                            "feature_frequency_per_", winName, "_chromosomeProfiles_",
                            "quantile", k, "_of_", quantiles,
                             "_by_log2_", libName, "_control_in_",
                             region, "_of_",
                             substring(featureName[1][1], first = 1, last = 5), "_in_",
                             paste0(substring(featureName, first = 10, last = 16),
                                    collapse = "_"), "_",
                             substring(featureName[1][1], first = 18), ".tsv"),
               row.names = F, col.names = T, quote = F, sep = "\t")
}
