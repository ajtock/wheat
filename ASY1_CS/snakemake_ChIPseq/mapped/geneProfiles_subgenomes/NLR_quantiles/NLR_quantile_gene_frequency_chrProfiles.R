#!/applications/R/R-3.5.0/bin/Rscript

# Calculate windowed feature frequencies for features within
# libName quantiles

# Usage:
# Rscript ./NLR_quantile_gene_frequency_chrProfiles.R cluster_members 'NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide' genes 4 1Mb 1000000 15

#libName <- "cluster_members"
#featureName <- unlist(strsplit("NLR_genes_in_Agenome_genomewide,NLR_genes_in_Bgenome_genomewide,NLR_genes_in_Dgenome_genomewide",
#                               split = ","))
#region <- "genes"
#quantiles <- 4
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15

args <- commandArgs(trailingOnly = T)
libName <- args[1]
featureName <- unlist(strsplit(args[2],
                               split = ","))
region <- args[3]
quantiles <- as.numeric(args[4])
winName <- args[5]
winSize <- as.numeric(args[6])
N <- as.numeric(args[7])

library(rtracklayer)
library(parallel)
library(GenomicRanges)

if(libName %in% c("cMMb", "cluster_members")) {
  outDir <- paste0("quantiles_by_", libName, "/")
} else {
  outDir <- paste0("quantiles_by_", sub("_\\w+", "", libName),
                   "_in_", region, "/")
}
plotDir <- paste0(outDir, "chrProfile_plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]

# Load table of features grouped into quantiles
if(libName %in% c("cMMb", "cluster_members")) {
  featuresDF <- read.table(paste0(outDir, "/WesternEurope/features_", quantiles, "quantiles",
                                  "_by_", libName, "_of_",
                                  substring(featureName[1][1], first = 1, last = 9), "_in_",
                                  paste0(substring(featureName, first = 14, last = 20),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 22), "_WesternEurope.txt"),
                           header = T, sep = "\t")
} else {
  featuresDF <- read.table(paste0(outDir, "/WesternEurope/features_", quantiles, "quantiles",
                                  "_by_", sub("_\\w+", "", libName), "_in_",
                                  region, "_of_",
                                  substring(featureName[1][1], first = 1, last = 9), "_in_",
                                  paste0(substring(featureName, first = 14, last = 20),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 22), "_WesternEurope.txt"),
                           header = T, sep = "\t")
}

quantilesDFlist <- mclapply(1:quantiles, function(k) {
  featuresDF[featuresDF$NLR_quantile == paste0("Quantile ", k),]
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

#library(doParallel)
#registerDoParallel(cores = quantiles)
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())

quantileProfiles <- NULL
for(k in 1:quantiles) {
  quantilekProfile <- NULL
  for(i in 1:length(chrs)) {
    # Count quantile features within adjacent windows
    featuresChr <- quantilesDFlist[[k]][quantilesDFlist[[k]]$seqnames == chrs[i],]
    windowsChrGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
    if(dim(featuresChr)[1] > 0) {
      featuresChrGR <- GRanges(seqnames = chrs[i],
                               ranges = IRanges(start = featuresChr$start,
                                                end = featuresChr$end),
                               strand = featuresChr$strand)
    } else {
      featuresChrGR <- GRanges()
    }
    winfeatures <- countOverlaps(query = windowsChrGR,
                                 subject = featuresChrGR,
                                 type = "any",
                                 ignore.strand = T)
    profileChr <- data.frame(chr = as.character(chrs[i]),
                             window = as.integer(start(windowsChrGR)),
                             features = as.integer(winfeatures),
                             stringsAsFactors = F)
    quantilekProfile <- rbind(quantilekProfile, profileChr)
  }
  quantileProfiles <- cbind(quantileProfiles, quantilekProfile[,3])
}
quantileProfiles <- data.frame(quantilekProfile[,1:2],
                               quantileProfiles,
                               stringsAsFactors = F)
colnames(quantileProfiles) <- c("chr", "window",
                                paste0("quantile", 1:quantiles))
if(libName %in% c("cMMb", "cluster_members")) {
  write.table(quantileProfiles,
              file = paste0(outDir,
                            "NLR_gene_frequency_per_", winName,
                            "_", quantiles, "quantiles_",
                             "_by_", libName, "_of_", 
                            substring(featureName[1][1], first = 1, last = 9), "_in_",
                            paste0(substring(featureName, first = 14, last = 20),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 22), "_unsmoothed.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
} else {
  write.table(quantileProfiles,
              file = paste0(outDir,
                            "NLR_gene_frequency_per_", winName,
                            "_", quantiles, "quantiles_",
                            "_by_log2_", libName, "_control_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 9), "_in_",
                            paste0(substring(featureName, first = 14, last = 20),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 22), "_unsmoothed.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
}
 
# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  quantileProfiles[quantileProfiles$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfiles <- mclapply(seq_along(chrProfiles), function(x) {
  filt_quantileProfiles <- NULL
  for(k in 1:quantiles) {
    filt_chrProfilek <- stats::filter(x = chrProfiles[[x]][,k+2],
                                      filter = f,
                                      sides = 2)
    filt_chrProfilek[1:flank] <- filt_chrProfilek[flank+1]
    filt_chrProfilek[(length(filt_chrProfilek)-flank+1):length(filt_chrProfilek)] <- filt_chrProfilek[(length(filt_chrProfilek)-flank)]
    filt_quantilekProfile <- data.frame(chr = as.character(chrProfiles[[x]]$chr),
                                        window = as.integer(chrProfiles[[x]]$window),
                                        filt_features = as.numeric(filt_chrProfilek),
                                        stringsAsFactors = F)
    filt_quantileProfiles <- cbind(filt_quantileProfiles, filt_quantilekProfile[,3])
  }
  filt_quantileProfiles <- data.frame(filt_quantilekProfile[,1:2],
                                      filt_quantileProfiles,
                                      stringsAsFactors = F)
  colnames(filt_quantileProfiles) <- c("chr", "window",
                                       paste0("filt_quantile", 1:quantiles))
  return(filt_quantileProfiles)
}, mc.cores = length(chrProfiles))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_quantileProfiles <- do.call(rbind, filt_chrProfiles)
if(libName %in% c("cMMb", "cluster_members")) {
  write.table(filt_quantileProfiles,
              file = paste0(outDir,
                            "NLR_gene_frequency_per_", winName,
                            "_", quantiles, "quantiles_",
                             "_by_", libName, "_of_", 
                            substring(featureName[1][1], first = 1, last = 9), "_in_",
                            paste0(substring(featureName, first = 14, last = 20),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 22), "_smoothed.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
} else {
  write.table(filt_quantileProfiles,
              file = paste0(outDir,
                            "NLR_gene_frequency_per_", winName,
                            "_", quantiles, "quantiles_",
                            "_by_log2_", libName, "_control_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 9), "_in_",
                            paste0(substring(featureName, first = 14, last = 20),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 22), "_smoothed.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
}
