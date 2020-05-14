#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 28.04.2020

# Create and plot correlation matrices of chromosome profiles for
# log2(ChIP/input), MNase, DNA methylation, and NLR-encoding and meiotic genes,
# genome-wide, subgenome-wide or separated into the following chromosome compartments a.k.a. partitions):
# 1. R1 and R3 (distal)
# 2. R2a and R2b (interstitial)
# 3. C (proximal)
# 4. heterochromatin (interstitial and proximal)
# 5. centromeres (defined by IWGSC (2018) Science 361 using CENH3 ChIP-seq data from Guo et al. (2016) PLOS Genet. 12)

# Usage:
# ./chrProfiles_correlation_matrices.R 1Mb 1000000 distal 'A,B,D'

#winName <- "1Mb"
#winSize <- 1000000
#region <- "distal"
#genomeName <- unlist(strsplit("A,B,D",
#                              split = ","))

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
region <- args[3]
genomeName <- unlist(strsplit(args[4],
                              split = ","))

library(GenomicRanges)
library(ggplot2)
library(ggcorrplot)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                               "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,2])
centromereEnd <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                             "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = chrStart,
                                     end = chrLens),
                    strand = "*")
if(length(genomeName) == 1) {
  genomeGR <- genomeGR[grep(genomeName,
                            seqnames(genomeGR))@values]
}

# Define region to be analysed
if(region == "distal") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "interstitial") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(chrPartitions$R1_R2a+1,
                                                 chrPartitions$C_R2b),
                                       end = c(chrPartitions$R2a_C,
                                               chrPartitions$R2b_R3-1)),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "proximal") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else {
  stop("region is not distal, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "distal") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "interstitial") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 3),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2a_C+1,
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrPartitions$C_R2b-1,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "proximal") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "genomewide") {
  maskGR <- GRanges()
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else {
  stop("region is not distal, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}


ASY1 <- read.table("log2_ASY1_CS_Rep1_ChIP_H3_input_SRR6350669_per_1Mb_smoothed.txt",
                   header = T)
ASY1GR <- GRanges(seqnames = ASY1$chr,
                  ranges = IRanges(start = ASY1$window,
                                   end = ASY1$window+winSize-1),
                  strand = "*",
                  value = ASY1$filt_log2ChIPcontrol)
ASY1GR_chrs <- lapply(seq_along(chrs), function(x) {
  ASY1GR_chr <- ASY1GR[seqnames(ASY1GR) == chrs[x]]
  end(ASY1GR_chr)[length(ASY1GR_chr)] <- chrLens[x]
  ASY1GR_chr
})
ASY1GR <- do.call(c, ASY1GR_chrs)

if(length(genomeName) == 1) {
  ASY1GR <- ASY1GR[grep(genomeName,
                        seqnames(ASY1GR))@values]
}
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_ASY1_overlap <- findOverlaps(query = maskGR,
                                  subject = ASY1GR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
ASY1GR <- ASY1GR[-subjectHits(mask_ASY1_overlap)]


#library(plyr)
#library(data.table)
#library(varhandle)
#library(zoo)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

## ChIPA profile
if(libNameChIPA %in% c("H3K4me3_ChIP_SRR6350668",
                       "H3K27me3_ChIP_SRR6350666",
                       "H3K36me3_ChIP_SRR6350670",
                       "H3K9ac_ChIP_SRR6350667",
                       "CENH3_ChIP_SRR1686799")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else if(libNameChIPA %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                              "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
}
profileChIPA <- read.table(paste0(covDirChIPA, libNameChIPA, "_MappedOn_wheat_v1.0_lowXM_",
                                  align, "_sort_norm_binSize", winName, ".bedgraph"))
if("chrUn" %in% chrs) {
  profileChIPA[dim(profileChIPA)[1],]$V3 <- chrLens[length(chrLens)]
} else {
  profileChIPA[dim(profileChIPA)[1],]$V3 <- 480980714
}
# Rows where the difference between end and start coordinates is > winSize
profileChIPA_bigWins <- profileChIPA[profileChIPA$V3-profileChIPA$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileChIPA <- profileChIPA[profileChIPA$V3-profileChIPA$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileChIPA_bigWinsList <- mclapply(seq_along(1:dim(profileChIPA_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileChIPA_bigWins[x,]$V2,
                      to = profileChIPA_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileChIPA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPA_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileChIPA_bigWins[x,]$V3)),
               V4 = as.numeric(profileChIPA_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileChIPA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPA_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileChIPA_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileChIPA_bigWinsDT <- rbindlist(profileChIPA_bigWinsList)
profileChIPA <- rbind.fill(profileChIPA, profileChIPA_bigWinsDT)
profileChIPA <- profileChIPA[order(profileChIPA$V1, profileChIPA$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPA <- profileChIPA[profileChIPA$V1 == chrs[x],]
  if(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V4))
  }
}, mc.cores = detectCores())
profileChIPA_chrLenValsADT <- rbindlist(chrLenValsAList)
profileChIPA <- rbind.fill(profileChIPA, profileChIPA_chrLenValsADT)
profileChIPA <- profileChIPA[order(profileChIPA$V1, profileChIPA$V2),]

profileChIPA <- data.frame(chr = as.character(profileChIPA$V1),
                           window = as.integer(profileChIPA$V2+1),
                           CPM = as.numeric(profileChIPA$V4),
                           stringsAsFactors = F)

## ControlA profile
if(libNameControlA == "MNase_Rep1") {
  covDirControlA <- paste0("/home/ajt200/analysis/wheat/",
                           "MNase/snakemake_ChIPseq/mapped/", align, "/bg/")
  profileControlA <- read.table(paste0(covDirControlA, "MNase_Rep1_MappedOn_wheat_v1.0_lowXM_",
                                       align, "_sort_norm_binSize", winName, ".bedgraph"))
} else if(libNameControlA == "H3_input_SRR6350669") {
  covDirControlA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                           "input/snakemake_ChIPseq/mapped/", align, "/bg/")
  profileControlA <- read.table(paste0(covDirControlA, "H3_input_SRR6350669_MappedOn_wheat_v1.0_lowXM_",
                                       align, "_sort_norm_binSize", winName, ".bedgraph"))
} else {
  if(!(libNameControlA %in% c("MNase_Rep1", "H3_input_SRR6350669"))) {
    stop("libNameControlA is neither MNase_Rep1 nor H3_input_SRR6350669")
  }
}
if("chrUn" %in% chrs) {
  profileControlA[dim(profileControlA)[1],]$V3 <- chrLens[length(chrLens)]
} else {
  profileControlA[dim(profileControlA)[1],]$V3 <- 480980714
}
# Rows where the difference between end and start coordinates is > winSize
profileControlA_bigWins <- profileControlA[profileControlA$V3-profileControlA$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileControlA <- profileControlA[profileControlA$V3-profileControlA$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
profileControlA_bigWinsList <- mclapply(seq_along(1:dim(profileControlA_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileControlA_bigWins[x,]$V2,
                      to = profileControlA_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileControlA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlA_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileControlA_bigWins[x,]$V3)),
               V4 = as.numeric(profileControlA_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileControlA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlA_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileControlA_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileControlA_bigWinsDT <- rbindlist(profileControlA_bigWinsList)
profileControlA <- rbind.fill(profileControlA, profileControlA_bigWinsDT)
profileControlA <- profileControlA[order(profileControlA$V1, profileControlA$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileControl <- profileControlA[profileControlA$V1 == chrs[x],]
  if(chrProfileControl[dim(chrProfileControl)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileControl[dim(chrProfileControl)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileControl[dim(chrProfileControl)[1],]$V4))
  }
}, mc.cores = detectCores())
profileControlA_chrLenValsADT <- rbindlist(chrLenValsAList)
profileControlA <- rbind.fill(profileControlA, profileControlA_chrLenValsADT)
profileControlA <- profileControlA[order(profileControlA$V1, profileControlA$V2),]

profileControlA <- data.frame(chr = as.character(profileControlA$V1),
                              window = as.integer(profileControlA$V2+1),
                              CPM = as.numeric(profileControlA$V4),
                              stringsAsFactors = F)

# Calculate log2((ChIP+1)/(Control+1)) coverage within each window
profileChIPAlog2 <- data.frame(chr = as.character(profileChIPA$chr),
                               window = as.integer(profileChIPA$window),
                               ChIP = as.numeric(profileChIPA$CPM),
                               control = as.numeric(profileControlA$CPM),
                               log2ChIPcontrol = as.numeric(log2((profileChIPA$CPM+1)/(profileControlA$CPM+1))),
                               stringsAsFactors = F)
                        
chrProfilesChIPA <- mclapply(seq_along(chrs), function(x) {
  profileChIPAlog2[profileChIPAlog2$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_chrProfilesChIPA <- mclapply(seq_along(chrProfilesChIPA), function(x) {
  # log2ChIPcontrol
  filt_log2ChIPcontrol <- stats::filter(x = chrProfilesChIPA[[x]]$log2ChIPcontrol,
                                        filter = f,
                                        sides = 2)
  filt_log2ChIPcontrol[1:flank] <- filt_log2ChIPcontrol[flank+1]
  filt_log2ChIPcontrol[(length(filt_log2ChIPcontrol)-flank+1):length(filt_log2ChIPcontrol)] <- filt_log2ChIPcontrol[(length(filt_log2ChIPcontrol)-flank)]

  # ChIP
  filt_ChIP <- stats::filter(x = chrProfilesChIPA[[x]]$ChIP,
                             filter = f,
                             sides = 2)
  filt_ChIP[1:flank] <- filt_ChIP[flank+1]
  filt_ChIP[(length(filt_ChIP)-flank+1):length(filt_ChIP)] <- filt_ChIP[(length(filt_ChIP)-flank)]

  # control
  filt_control <- stats::filter(x = chrProfilesChIPA[[x]]$control,
                                        filter = f,
                                        sides = 2)
  filt_control[1:flank] <- filt_control[flank+1]
  filt_control[(length(filt_control)-flank+1):length(filt_control)] <- filt_control[(length(filt_control)-flank)]

  data.frame(chr = as.character(chrProfilesChIPA[[x]]$chr),
             window = as.integer(chrProfilesChIPA[[x]]$window),
             filt_ChIP = as.numeric(filt_ChIP),
             filt_control = as.numeric(filt_control),
             filt_log2ChIPcontrol = as.numeric(filt_log2ChIPcontrol),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesChIPA))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_profileChIPAlog2 <- do.call(rbind, filt_chrProfilesChIPA)

write.table(profileChIPAlog2,
            file = paste0("log2_", libNameChIPA, "_", libNameControlA, "_per_", winName, "_unsmoothed.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(filt_profileChIPAlog2,
            file = paste0("log2_", libNameChIPA, "_", libNameControlA, "_per_", winName, "_smoothed.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

