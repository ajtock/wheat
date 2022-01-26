#!/applications/R/R-4.0.0/bin/Rscript

# Compute signals of epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq,
# within physical marker intervals to enable correlations with wild type and fancm mutant
# crossover rates (cM/Mb)

# Usage:
# ./AxC_markers-mapping-summary_mapped_markers.R ASY1_CS ASY1_CS_Rep1_ChIP input H3_input_SRR6350669 both

#markChIP <- "ASY1_CS"
#libNameChIP <- "ASY1_CS_Rep1_ChIP"
#markControl <- "input"
#libNameControl <- "H3_input_SRR6350669"
#align <- both

args <- commandArgs(trailingOnly = T)
markChIP <- args[1]
libNameChIP <- args[2]
markControl <- args[3]
libNameControl <- args[4]
align <- args[5]

options(stringsAsFactors = F)
library(data.table)
library(parallel)
library(GenomicRanges)

# Genomic definitions
winSize <- 1
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
                                           end = c(tab_chr$pos)),
                          strand = "*",
                          physical_marker = tab_chr$physical_marker,
                          wt_marker = tab_chr$control_marker,
                          wt_cM = tab_chr$control_cM,
                          fancm_marker = tab_chr$fancm_marker,
                          fancm_cM = tab_chr$fancm_cM) 
  interGR <- c(interGR, inter_chr_GR)
}

## ChIP profile
if(libNameChIP %in% c("H3K4me3_ChIP_SRR6350668",
                     "H3K27me3_ChIP_SRR6350666",
                     "H3K36me3_ChIP_SRR6350670",
                     "H3K9ac_ChIP_SRR6350667",
                     "CENH3_ChIP_SRR1686799",
                     "H3_input_SRR6350669")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                       markChIP, "/snakemake_ChIPseq_MAPQ0_XM6_normNone/mapped/",
                       align, "/bg/")
} else if(libNameChIP %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                             "H3K27ac_Rep1_ChIP_SRR8126621")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
                       markChIP, "/snakemake_ChIPseq_MAPQ0_XM6_normNone/mapped/",
                       align, "/bg/")
} else {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/",
                       markChIP, "/snakemake_ChIPseq_MAPQ0_XM6_normNone/mapped/",
                       align, "/bg/")
}
ChIP <- read.table(paste0(covDirChIP, libNameChIP, "_MappedOn_wheat_v1.0_lowXM_",
                          align, "_sort_norm.perbase"))

# Load table of normalised ChIP coverage values             
#ASY1 <- read.table("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/bg/ASY1_CS_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm.perbase",
#                   header = F)
#ChIP <- ASY1

# Function to convert coverage table into GRanges
makeGR <- function(bedgraph) {
  GRanges(seqnames = bedgraph[,1],
          ranges = IRanges(start = bedgraph[,2],
                           end = bedgraph[,2]),
          strand = "*",
          val = bedgraph[,3])
}

ChIP_GR <- makeGR(bedgraph = ChIP)

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

inter_ChIP_GR <- ChIP_GR[subjectHits(fOverlaps_ChIP)]
#inter_ChIP_GR_win <- inter_ChIP_GR[width(inter_ChIP_GR) > 1]

#inter_ChIP_GR_exp_list <- mclapply(1:length(inter_ChIP_GR), function(x) {
inter_ChIP_GR_exp_list <- lapply(1:length(inter_ChIP_GR), function(x) {
  print(x)
  tmp <- rep(inter_ChIP_GR[x], width(inter_ChIP_GR[x]))
  for(i in 1:width(inter_ChIP_GR[x])) {
    ranges(tmp[i]) <- IRanges(start = start(tmp[i]) + i - 1,
                              end = start(tmp[i]) + i - 1)
  }
  tmp
})
#}, mc.cores = detectCores(), mc.preschedule = T)




## Control profile
if(libNameControl == "MNase_Rep1") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/",
                          "MNase/snakemake_ChIPseq/mapped/", align, "/bg/")
  Control <- read.table(paste0(covDirControl, "MNase_Rep1_MappedOn_wheat_v1.0_lowXM_",
                               align, "_sort_norm_binSize", winName, ".perbase"))
} else if(libNameControl == "H3_input_SRR6350669") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                          "input/snakemake_ChIPseq/mapped/", align, "/bg/")
  Control <- read.table(paste0(covDirControl, "H3_input_SRR6350669_MappedOn_wheat_v1.0_lowXM_",
                               align, "_sort_norm_binSize", winName, ".perbase"))
} else {
  if(!(libNameControl %in% c("MNase_Rep1", "H3_input_SRR6350669"))) {
    stop("libNameControl is neither MNase_Rep1 nor H3_input_SRR6350669")
  }
}


chr_tabGR_str_x <- chr_tabGR_str[subjectHits(fOverlaps_str[queryHits(fOverlaps_str) == featNum])]



                   

# Rows where the difference between end and start coordinates is > winSize
ASY1_bigWins <- ASY1[ASY1$V3-ASY1$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
ASY1 <- ASY1[ASY1$V3-ASY1$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end
ASY1_bigWinsList <- mclapply(seq_along(1:dim(ASY1_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = ASY1_bigWins[x,]$V2,
                      to = ASY1_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < ASY1_bigWins[x,]$V3) {
    data.frame(V1 = as.character(ASY1_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 ASY1_bigWins[x,]$V3)),
               V4 = as.numeric(ASY1_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == ASY1_bigWins[x,]$V3) {
    data.frame(V1 = as.character(ASY1_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(ASY1_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())
 
ASY1_bigWinsDT <- rbindlist(ASY1_bigWinsList)
ASY1 <- rbind.fill(ASY1, ASY1_bigWinsDT)
ASY1 <- ASY1[order(ASY1$V1, ASY1$V2),]

chrLenValsAList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPA <- ASY1[ASY1$V1 == chrs[x],]
  if(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V4))
  }
}, mc.cores = detectCores())
ASY1_chrLenValsADT <- rbindlist(chrLenValsAList)
ASY1 <- rbind.fill(ASY1, ASY1_chrLenValsADT)
ASY1 <- ASY1[order(ASY1$V1, ASY1$V2),]

ASY1 <- data.frame(chr = as.character(ASY1$V1),
                           window = as.integer(ASY1$V2+1),
                           CPM = as.numeric(ASY1$V4),
                           stringsAsFactors = F)


