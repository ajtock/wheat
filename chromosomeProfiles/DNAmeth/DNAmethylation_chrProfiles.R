#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# csmit -m 300G -c 47 "Rscript ./DNAmethylation_chrProfiles.R BSseq_Rep8a_SRR6792678 1Mb 1000000 15"

#DNAmethName <- "BSseq_Rep8a_SRR6792678"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15

args <- commandArgs(trailingOnly = T)
DNAmethName <- args[1]
winName <- args[2]
winSize <- as.numeric(args[3])
N <- as.numeric(args[4])

library(parallel)
library(GenomicRanges)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
#chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
#chrLens <- chrLens[-length(chrLens)]

if(grepl(pattern = "SRR67926", x = DNAmethName)) {
  DNAmethDir <- "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/BSseq/snakemake_BSseq/coverage/"
} else {
  stop(paste0("DNAmethName is not compatible with the specified coverage path"))
}

# Load processed DNA methylation data
mCG <- read.table(gzfile(paste0(DNAmethDir,
                                DNAmethName, "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_",
                                "CpG.gz.bismark.cov.gz")),
                  header = F, stringsAsFactors = F,
                  colClasses = c(rep(NA, 2), "NULL", NA, rep("NULL", 2)))
colnames(mCG) <- c("chr", "pos", "propMeth")
mCG$propMeth <- as.numeric(mCG$propMeth/100)
mCHG <- read.table(gzfile(paste0(DNAmethDir,
                                 DNAmethName, "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_",
                                 "CHG.gz.bismark.cov.gz")),
                   header = F, stringsAsFactors = F,
                   colClasses = c(rep(NA, 2), "NULL", NA, rep("NULL", 2)))
colnames(mCHG) <- c("chr", "pos", "propMeth")
mCHG$propMeth <- as.numeric(mCHG$propMeth/100)
mCHH <- read.table(gzfile(paste0(DNAmethDir,
                                 DNAmethName, "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_",
                                 "CHH.gz.bismark.cov.gz")),
                   header = F, stringsAsFactors = F,
                   colClasses = c(rep(NA, 2), "NULL", NA, rep("NULL", 2)))
colnames(mCHH) <- c("chr", "pos", "propMeth")
mCHH$propMeth <- as.numeric(mCHH$propMeth/100)

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

profileDNAmeth <- NULL
for(i in 1:length(chrs)) {
  print(chrs[i])
  chr_windowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  # mCG
  chr_mCG <- mCG[mCG$chr == chrs[i],]
  chr_mCG_GR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chr_mCG$pos,
                                         width = 1),
                        strand = "*")
  # Identify overlapping windows and cytosine coordinates in CG context
  mCG_fOverlaps <- findOverlaps(query = chr_windowsGR,
                                subject = chr_mCG_GR,
                                type = "any",
                                select = "all",
                                ignore.strand = TRUE)
  # Convert mCG_fOverlaps into list object equivalent to that
  # generated by segmentSeq::getOverlaps(), in which each
  # list element corresponds to a genomic window of
  # sequentially numbered cytosine coordinates
  mCG_fOverlapsList <- mclapply(seq_along(unique(queryHits(mCG_fOverlaps))),
                       function(x) {
                         subjectHits(mCG_fOverlaps[queryHits(mCG_fOverlaps) == x])
                       }, mc.cores = detectCores(), mc.preschedule = T)  
  # Calculate mean methylation proportion in each window
  win_mCG <- sapply(mCG_fOverlapsList, function(x) {
               mean(chr_mCG$propMeth[x])
             })

  # mCHG
  chr_mCHG <- mCHG[mCHG$chr == chrs[i],]
  chr_mCHG_GR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = chr_mCHG$pos,
                                          width = 1),
                         strand = "*")
  # Identify overlapping windows and cytosine coordinates in CHG context
  mCHG_fOverlaps <- findOverlaps(query = chr_windowsGR,
                                 subject = chr_mCHG_GR,
                                 type = "any",
                                 select = "all",
                                 ignore.strand = TRUE)
  # Convert mCHG_fOverlaps into list object equivalent to that
  # generated by segmentSeq::getOverlaps(), in which each
  # list element corresponds to a genomic window of
  # sequentially numbered cytosine coordinates
  mCHG_fOverlapsList <- mclapply(seq_along(unique(queryHits(mCHG_fOverlaps))),
                        function(x) {
                          subjectHits(mCHG_fOverlaps[queryHits(mCHG_fOverlaps) == x])
                        }, mc.cores = detectCores(), mc.preschedule = T)  
  # Calculate mean methylation proportion in each window
  win_mCHG <- sapply(mCHG_fOverlapsList, function(x) {
                mean(chr_mCHG$propMeth[x])
              })

  # mCHH
  chr_mCHH <- mCHH[mCHH$chr == chrs[i],]
  chr_mCHH_GR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = chr_mCHH$pos,
                                          width = 1),
                         strand = "*")
  # Identify overlapping windows and cytosine coordinates in CHH context
  mCHH_fOverlaps <- findOverlaps(query = chr_windowsGR,
                                 subject = chr_mCHH_GR,
                                 type = "any",
                                 select = "all",
                                 ignore.strand = TRUE)
  # Convert mCHH_fOverlaps into list object equivalent to that
  # generated by segmentSeq::getOverlaps(), in which each
  # list element corresponds to a genomic window of
  # sequentially numbered cytosine coordinates
  mCHH_fOverlapsList <- mclapply(seq_along(unique(queryHits(mCHH_fOverlaps))),
                        function(x) {
                          subjectHits(mCHH_fOverlaps[queryHits(mCHH_fOverlaps) == x])
                        }, mc.cores = detectCores(), mc.preschedule = T)  
  # Calculate mean methylation proportion in each window
  win_mCHH <- sapply(mCHH_fOverlapsList, function(x) {
                mean(chr_mCHH$propMeth[x])
              })

  # Mean of all 3 contexts
  win_mC <- sapply(seq_along(win_mCG), function(x) {
              mean(c(win_mCG[x], win_mCHG[x], win_mCHH[x]))
            })

  # Make data.frame
  chr_profile <- data.frame(chr = chrs[i],
                            window = start(chr_windowsGR),
                            mCG = win_mCG,
                            mCHG = win_mCHG,
                            mCHH = win_mCHH,
                            mC = win_mC,
                            stringsAsFactors = F)
  # Iteratively combine profiles from each chromosome
  # in a data.frame (profileDNAmeth) that grows with the
  # completion of each loop
  profileDNAmeth <- rbind(profileDNAmeth, chr_profile)
} 
write.table(profileDNAmeth,
            file = paste0(DNAmethName, "_per_", winName, ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

chrProfilesDNAmeth <- mclapply(seq_along(chrs) , function(x) {
  profileDNAmeth[profileDNAmeth$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfilesDNAmeth <- mclapply(seq_along(chrProfilesDNAmeth), function(x) {
  # mCG
  filt_chrProfile_mCG <- stats::filter(x = chrProfilesDNAmeth[[x]]$mCG,
                                       filter = f,
                                       sides = 2)
  filt_chrProfile_mCG[1:flank] <- filt_chrProfile_mCG[flank+1]
  filt_chrProfile_mCG[(length(filt_chrProfile_mCG)-flank+1):length(filt_chrProfile_mCG)] <- filt_chrProfile_mCG[(length(filt_chrProfile_mCG)-flank)]
  # mCHG
  filt_chrProfile_mCHG <- stats::filter(x = chrProfilesDNAmeth[[x]]$mCHG,
                                        filter = f,
                                        sides = 2)
  filt_chrProfile_mCHG[1:flank] <- filt_chrProfile_mCHG[flank+1]
  filt_chrProfile_mCHG[(length(filt_chrProfile_mCHG)-flank+1):length(filt_chrProfile_mCHG)] <- filt_chrProfile_mCHG[(length(filt_chrProfile_mCHG)-flank)]
  # mCHH
  filt_chrProfile_mCHH <- stats::filter(x = chrProfilesDNAmeth[[x]]$mCHH,
                                        filter = f,
                                        sides = 2)
  filt_chrProfile_mCHH[1:flank] <- filt_chrProfile_mCHH[flank+1]
  filt_chrProfile_mCHH[(length(filt_chrProfile_mCHH)-flank+1):length(filt_chrProfile_mCHH)] <- filt_chrProfile_mCHH[(length(filt_chrProfile_mCHH)-flank)]
  # mC
  filt_chrProfile_mC <- stats::filter(x = chrProfilesDNAmeth[[x]]$mC,
                                      filter = f,
                                      sides = 2)
  filt_chrProfile_mC[1:flank] <- filt_chrProfile_mC[flank+1]
  filt_chrProfile_mC[(length(filt_chrProfile_mC)-flank+1):length(filt_chrProfile_mC)] <- filt_chrProfile_mC[(length(filt_chrProfile_mC)-flank)]
  # Combine in data.frame
  data.frame(chr = as.character(chrProfilesDNAmeth[[x]]$chr),
             window = as.integer(chrProfilesDNAmeth[[x]]$window),
             filt_mCG = as.numeric(filt_chrProfile_mCG),
             filt_mCHG = as.numeric(filt_chrProfile_mCHG),
             filt_mCHH = as.numeric(filt_chrProfile_mCHH),
             filt_mC = as.numeric(filt_chrProfile_mC),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesDNAmeth))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_profileDNAmeth <- do.call(rbind, filt_chrProfilesDNAmeth)

write.table(filt_profileDNAmeth,
            file = paste0(DNAmethName, "_per_", winName, "_smoothed.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
