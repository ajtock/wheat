#!/applications/R/R-4.0.0/bin/Rscript

# Compute signals of epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq,
# within physical marker intervals to enable correlations with wild type and fancm mutant
# crossover rates (cM/Mb)

# Usage:
# ./append_DNAmeth.R BSseq_Rep8a_SRR6792678 CpG

#DNAmethName <- "BSseq_Rep8a_SRR6792678"
#context <- "CpG"

args <- commandArgs(trailingOnly = T)
DNAmethName <- args[1]
context <- args[2] 

options(stringsAsFactors = F)
library(data.table)
library(parallel)
library(GenomicRanges)
library(dplyr)
library(plyr)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

# Table of inter-marker cM values in WT and fancm
#tab <- fread("AxC_markers-mapping-summary_mapped_markers.tsv")
tab <- read.table("AxC_markers-mapping-summary_mapped_markers.tsv",
                  header = T)
tab$pos <- as.integer(gsub(",", "", tab$pos))
tab$chr <- sub("^", "chr", tab$chr)
tab <- tab[ with(tab, base::order(chr, pos)) , ]

interGR <- GRanges()
for(i in unique(tab$chr)) {
  print(i)
  tab_chr <- tab[tab$chr == i,]
  inter_chr_GR <- GRanges(seqnames = c(tab_chr$chr),
                          ranges = IRanges(start = c(1, tab_chr[-nrow(tab_chr),]$pos+1),
                                           end = tab_chr$pos),
                          strand = "*",
                          physical_marker = tab_chr$physical_marker,
                          wt_marker = tab_chr$control_marker,
                          wt_cM = tab_chr$control_cM,
                          fancm_marker = tab_chr$fancm_marker,
                          fancm_cM = tab_chr$fancm_cM) 
  interGR <- c(interGR, inter_chr_GR)
}

#interGR_poswidths <- interGR[width(interGR) > 0]
interGR[width(interGR) == 0]
quantile(width(interGR), probs = seq(0, 1, 0.05))
#          0%           5%          10%          15%          20%          25% 
#        0.00        86.80       263.60       773.20      2502.60      7338.50 
#         30%          35%          40%          45%          50%          55% 
#    22745.90     61824.75    101485.60    145696.75    214491.50    297267.35 
#         60%          65%          70%          75%          80%          85% 
#   413309.40    532912.60    703690.50    970857.50   1340466.00   1869535.45 
#         90%          95%         100% 
#  3301806.00   8655394.25 703285697.00 
#pdf("marker_interval_widths_hist.pdf")
#hist(width(interGR), breaks = 100000, xlim = c(0, 1e6))
#dev.off()

if(grepl(pattern = "SRR67926", x = DNAmethName)) {
  DNAmethDir <- "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/BSseq/snakemake_BSseq/coverage/"
} else {
  stop(paste0("DNAmethName is not compatible with the specified coverage path"))
}

# Load processed DNA methylation data
mC <- read.table(gzfile(paste0(DNAmethDir,
                                DNAmethName, "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_",
                                context, ".gz.bismark.cov.gz")),
                  header = F, stringsAsFactors = F,
                  colClasses = c(rep(NA, 2), "NULL", NA, rep("NULL", 2)))
colnames(mC) <- c("chr", "pos", "propMeth")
mC$propMeth <- as.numeric(mC$propMeth/100)


# Function to convert coverage table into GRanges
makeGR <- function(bedgraph) {
  GRanges(seqnames = bedgraph[,1],
          ranges = IRanges(start = bedgraph[,2],
                           end = bedgraph[,2]),
          strand = "*",
          val = bedgraph[,3])
}

mC_GR <- makeGR(bedgraph = mC)

# Function to find coverage windows that overlap marker intervals
fOverlaps <- function(interGR, datGR) {
  fOverlaps_obj <- findOverlaps(query = interGR,
                                subject = datGR,
                                type = "any",
                                select = "all",
                                ignore.strand = T)
  fOverlaps_obj
}

#fOverlaps_mC <- fOverlaps(interGR = interGR, datGR = mC_GR)

profileDNAmeth <- NULL
for(chrName in seqlevels(interGR)) {
  print(chrName)
  chr_interGR <- interGR[seqnames(interGR) == chrName]

  # mC
  chr_mC_GR <- mC_GR[seqnames(mC_GR) == chrName]

  # Identify overlapping windows and cytosine coordinates in CG context
  fOverlaps_mC <- fOverlaps(interGR = chr_interGR, datGR = chr_mC_GR)

  # Convert fOverlaps_mC* into list object equivalent to that
  # generated by segmentSeq::getOverlaps(), in which each
  # list element corresponds to a genomic window (marker interval) of
  # sequentially numbered overlapping cytosine coordinates
  fOverlaps_mC_list <- mclapply(seq_along(unique(queryHits(fOverlaps_mC))),
                          function(x) {
                            subjectHits(fOverlaps_mC[queryHits(fOverlaps_mC) == x])
                          }, mc.cores = detectCores(), mc.preschedule = T)
#  mC_fOverlapsList <- mclapply(seq_along(unique(queryHits(mC_fOverlaps))),
#                         function(x) {
#                           subjectHits(mC_fOverlaps[queryHits(mC_fOverlaps) == x])
#                         }, mc.cores = detectCores(), mc.preschedule = T)
  # Calculate mean methylation proportion in each window
  win_mC <- sapply(mC_fOverlapsList, function(x) {
               mean(chr_mC$propMeth[x])
             })


# Convert fOverlaps_mC* into list object equivalent to that
# generated by segmentSeq::getOverlaps(), in which each
# list element corresponds to a genomic window (marker interval) of
# sequentially numbered overlapping cytosine coordinates

library(segmentSeq)
fOverlaps_mC_list <- segmentSeq::getOverlaps(coordinates = interGR,
                                              segments = mC_GR,
                                              overlapType = "overlapping",
                                              whichOverlaps = T,
                                              ignoreStrand = T)

#fOverlaps_mC_list <- mclapply(seq_along(unique(queryHits(fOverlaps_mC))),
fOverlaps_mC_list <- lapply(seq_along(unique(queryHits(fOverlaps_mC))),
                        function(x) {
                          print(x)
                          subjectHits(fOverlaps_mC[queryHits(fOverlaps_mC) == x])
                        })
#                        }, mc.cores = detectCores(), mc.preschedule = T)
fOverlaps_mCHG_list <- mclapply(seq_along(unique(queryHits(fOverlaps_mCHG))),
                         function(x) {
                           subjectHits(fOverlaps_mCHG[queryHits(fOverlaps_mCHG) == x])
                         }, mc.cores = detectCores(), mc.preschedule = T)
# Calculate mean methylation proportion in each window
win_mCHG <- sapply(fOverlaps_mCHG_list, function(x) {
              mean(chr_mCHG$propMeth[x], na.rm = T)
            })






## ChIP profile
if(libNameChIP %in% c("H3K4me3_ChIP_SRR6350668",
                      "H3K27me3_ChIP_SRR6350666",
                      "H3K36me3_ChIP_SRR6350670",
                      "H3K9ac_ChIP_SRR6350667",
                      "CENH3_ChIP_SRR1686799",
                      "input_SRR6350669")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                       markChIP, "/snakemake_ChIPseq/mapped/",
                       align, "/bg/")
} else if(libNameChIP %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                             "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                       markChIP, "/snakemake_ChIPseq/mapped/",
                       align, "/bg/")
} else {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/",
                       markChIP, "/snakemake_ChIPseq/mapped/",
                       align, "/bg/")
}
ChIP <- read.table(paste0(covDirChIP, libNameChIP, "_MappedOn_wheat_v1.0_lowXM_",
                          align, "_sort_norm_binSize", genomeBinName , ".bedgraph"))

# Rows where the difference between end and start coordinates is > genomeBinSize
ChIP_bigWins <- ChIP[ChIP$V3-ChIP$V2 > genomeBinSize,]
# Rows where the difference between end and start coordinates is == genomeBinSize
ChIP <- ChIP[ChIP$V3-ChIP$V2 == genomeBinSize,]

# Create a list of big windows, each split into windows of genomeBinSize,
# or < genomeBinSize if at chromosome end
ChIP_bigWinsList <- mclapply(seq_along(1:dim(ChIP_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = ChIP_bigWins[x,]$V2,
                      to = ChIP_bigWins[x,]$V3,
                      by = genomeBinSize)

  if(bigWinsSplit[length(bigWinsSplit)] < ChIP_bigWins[x,]$V3) {
    data.frame(V1 = as.character(ChIP_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize,
                                 ChIP_bigWins[x,]$V3)),
               V4 = as.numeric(ChIP_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == ChIP_bigWins[x,]$V3) {
    data.frame(V1 = as.character(ChIP_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize),
               V4 = as.numeric(ChIP_bigWins[x,]$V4))
  }
}, mc.cores = detectCores(), mc.preschedule = T)

ChIP_bigWinsDT <- rbindlist(ChIP_bigWinsList)
ChIP <- rbind.fill(ChIP, ChIP_bigWinsDT)
ChIP <- ChIP[order(ChIP$V1, ChIP$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIP <- ChIP[ChIP$V1 == chrs[x],]
  if(chrProfileChIP[dim(chrProfileChIP)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIP[dim(chrProfileChIP)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIP[dim(chrProfileChIP)[1],]$V4))
  }
}, mc.cores = detectCores(), mc.preschedule = F)
ChIP_chrLenValsDT <- rbindlist(chrLenValsList)
ChIP <- rbind.fill(ChIP, ChIP_chrLenValsDT)
ChIP <- ChIP[order(ChIP$V1, ChIP$V2),]

colnames(ChIP) <- c("chr", "start", "end", "val")

ChIP$start <- ChIP$start+1


## Control profile
if(libNameControl == "MNase_Rep1") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/",
                          "MNase/snakemake_ChIPseq/mapped/", align, "/bg/")
  Control <- read.table(paste0(covDirControl, "MNase_Rep1_MappedOn_wheat_v1.0_lowXM_",
                               align, "_sort_norm_binSize", genomeBinName, ".bedgraph"))
} else if(libNameControl == "input_SRR6350669") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                          "input/snakemake_ChIPseq/mapped/", align, "/bg/")
  Control <- read.table(paste0(covDirControl, "input_SRR6350669_MappedOn_wheat_v1.0_lowXM_",
                               align, "_sort_norm_binSize", genomeBinName, ".bedgraph"))
} else {
  if(!(libNameControl %in% c("MNase_Rep1", "input_SRR6350669"))) {
    stop("libNameControl is neither MNase_Rep1 nor input_SRR6350669")
  }
}
Control <- read.table(paste0(covDirControl, libNameControl, "_MappedOn_wheat_v1.0_lowXM_",
                             align, "_sort_norm_binSize", genomeBinName , ".bedgraph"))

# Rows where the difference between end and start coordinates is > genomeBinSize
Control_bigWins <- Control[Control$V3-Control$V2 > genomeBinSize,]
# Rows where the difference between end and start coordinates is == genomeBinSize
Control <- Control[Control$V3-Control$V2 == genomeBinSize,]

# Create a list of big windows, each split into windows of genomeBinSize,
# or < genomeBinSize if at chromosome end
Control_bigWinsList <- mclapply(seq_along(1:dim(Control_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = Control_bigWins[x,]$V2,
                      to = Control_bigWins[x,]$V3,
                      by = genomeBinSize)

  if(bigWinsSplit[length(bigWinsSplit)] < Control_bigWins[x,]$V3) {
    data.frame(V1 = as.character(Control_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize,
                                 Control_bigWins[x,]$V3)),
               V4 = as.numeric(Control_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == Control_bigWins[x,]$V3) {
    data.frame(V1 = as.character(Control_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize),
               V4 = as.numeric(Control_bigWins[x,]$V4))
  }
}, mc.cores = detectCores(), mc.preschedule = T)

Control_bigWinsDT <- rbindlist(Control_bigWinsList)
Control <- rbind.fill(Control, Control_bigWinsDT)
Control <- Control[order(Control$V1, Control$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileControl <- Control[Control$V1 == chrs[x],]
  if(chrProfileControl[dim(chrProfileControl)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileControl[dim(chrProfileControl)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileControl[dim(chrProfileControl)[1],]$V4))
  }
}, mc.cores = detectCores(), mc.preschedule = F)
Control_chrLenValsDT <- rbindlist(chrLenValsList)
Control <- rbind.fill(Control, Control_chrLenValsDT)
Control <- Control[order(Control$V1, Control$V2),]

colnames(Control) <- c("chr", "start", "end", "val")

Control$start <- Control$start+1

stopifnot(identical(ChIP$chr, Control$chr))
stopifnot(identical(ChIP$start, Control$start))
stopifnot(identical(ChIP$end, Control$end))


# Function to convert coverage table into GRanges
makeGR <- function(bedgraph) {
  GRanges(seqnames = bedgraph[,1],
          ranges = IRanges(start = bedgraph[,2],
                           end = bedgraph[,3]),
          strand = "*",
          val = bedgraph[,4])
}

ChIP_GR <- makeGR(bedgraph = ChIP)

Control_GR <- makeGR(bedgraph = Control)

# Function to find coverage windows that overlap marker intervals
fOverlaps <- function(interGR, datGR) {
  fOverlaps_obj <- findOverlaps(query = interGR,
                                subject = datGR,
                                type = "any",
                                select = "all",
                                ignore.strand = T)
  fOverlaps_obj
}

fOverlaps_ChIP <- fOverlaps(interGR = interGR, datGR = ChIP_GR)

fOverlaps_Control <- fOverlaps(interGR = interGR, datGR = Control_GR)

# Function to calculate average per-base values for a given interval x
makeDF_x <- function(fOverlaps_ChIP, fOverlaps_Control, ChIP_GR, Control_GR, interGR, interNum) {

  ChIP_GR_x <- ChIP_GR[subjectHits(fOverlaps_ChIP[queryHits(fOverlaps_ChIP) == interNum])]
#  ChIP_GR_x_VPB <- sum(ChIP_GR_x$val) / (sum(width(ChIP_GR_x))/1e3)
  ChIP_GR_x_VPB <- mean(ChIP_GR_x$val, na.rm = T) 

  Control_GR_x <- Control_GR[subjectHits(fOverlaps_Control[queryHits(fOverlaps_Control) == interNum])]
#  Control_GR_x_VPB <- sum(Control_GR_x$val) / (sum(width(Control_GR_x))/1e3)
  Control_GR_x_VPB <- mean(Control_GR_x$val, na.rm = T) 

  data.frame(interGR[interNum],
             ChIP = ChIP_GR_x_VPB,
             Control = Control_GR_x_VPB,
             log2val = log2( (ChIP_GR_x_VPB+1) / (Control_GR_x_VPB+1) ) )
}

# Apply makeDF_x function to each range in interGR
makeDF_x_list <- mclapply(1:length(interGR), function(x) {
  makeDF_x(fOverlaps_ChIP = fOverlaps_ChIP,
           fOverlaps_Control = fOverlaps_Control,
           ChIP_GR = ChIP_GR,
           Control_GR = Control_GR,
           interGR = interGR,
           interNum = x)
}, mc.cores = detectCores(), mc.preschedule = T)

makeDF <- dplyr::bind_rows(makeDF_x_list, .id = "column_label")

makeDF <- makeDF[,c(-1, -6)]

colnames(makeDF) <- c("chr",
                      colnames(makeDF)[2:9],
                      markChIP, markControl,
                      paste0("log2_", markChIP, "_", markControl))

write.table(makeDF,
            file = paste0("AxC_mapped_marker_intervals_",
                          libNameChIP, "_", libNameControl, "_",
                          align, "_binSize", genomeBinName, ".tsv"),
             quote = F, sep = "\t", row.names = F, col.names = T)


##inter_ChIP_GR <- ChIP_GR[subjectHits(fOverlaps_ChIP)]
###inter_ChIP_GR_win <- inter_ChIP_GR[width(inter_ChIP_GR) > 1]
##
###inter_ChIP_GR_exp_list <- mclapply(1:length(inter_ChIP_GR), function(x) {
##inter_ChIP_GR_exp_list <- lapply(1:length(inter_ChIP_GR), function(x) {
##  print(x)
##  tmp <- rep(inter_ChIP_GR[x], width(inter_ChIP_GR[x]))
##  for(i in 1:width(inter_ChIP_GR[x])) {
##    ranges(tmp[i]) <- IRanges(start = start(tmp[i]) + i - 1,
##                              end = start(tmp[i]) + i - 1)
##  }
##  tmp
##})
###}, mc.cores = detectCores(), mc.preschedule = T)
