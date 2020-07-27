#!/applications/R/R-3.5.0/bin/Rscript

# Plot widths of peaks and ranLoc that overlap different gene regions
# (bodies, promoters, terminators, TSS+500 bp, TTS-500 bp)

# Usage:
# ./DMC1_Rep1_ChIP_peaks_ranLoc_widths_overlapping_gene_regions.R DMC1_Rep1_ChIP DMC1 'chr3B'

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
library(ggplot2)

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
                    R2C[x], "_randomLoci.gff"),
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


# Load TEs
otherNames <- c(
                "CACTA_DTC",
                "Harbinger_DTH",
                "hAT_DTA",
                "Helitrons_DHH",
                "Mariner_DTT",
                "Mutator_DTM",
                "MITE_DMI",
                "Unclassified_class_2_DXX",
                "Unclassified_with_TIRs_DTX",
                "Copia_LTR_RLC",
                "Gypsy_LTR_RLG",
                "LINE_RIX",
                "SINE_SIX",
                "Unclassified_LTR_RLX",
                "Unclassified_repeats_XXX"
               )

TEdir <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/superfamilies/"

othersGRL <- lapply(seq_along(otherNames), function(x) {
  TEs <- read.table(paste0(TEdir,
                           "iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
                           otherNames[x], ".bed"), header = F)
  TEsGRtmp <- GRanges(seqnames = TEs$V1,
                      ranges = IRanges(start = TEs$V2+1,
                                       end = TEs$V3),
                      strand = "*")
  TEsGRtmp
})

CACTA_DTC_GR <- othersGRL[[1]]
Copia_LTR_RLC_GR <- othersGRL[[10]]
Gypsy_LTR_RLG_GR <- othersGRL[[11]]

# Load introns
intronsA <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_introns_in_Agenome_genomewide.gff3",
                    colClasses = c(NA,
                                   rep("NULL", 2),
                                   rep(NA, 2),
                                   "NULL", NA, "NULL", NA, "NULL"))
intronsB <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_introns_in_Bgenome_genomewide.gff3",
                    colClasses = c(NA,
                                   rep("NULL", 2),
                                   rep(NA, 2),
                                   "NULL", NA, "NULL", NA, "NULL"))
intronsD <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_introns_in_Dgenome_genomewide.gff3",
                    colClasses = c(NA,
                                   rep("NULL", 2),
                                   rep(NA, 2),
                                   "NULL", NA, "NULL", NA, "NULL"))
introns <- do.call(rbind, list(intronsA, intronsB, intronsD))

colnames(introns) <- c("chr", "start", "end", "strand", "intronID")
introns <- introns[introns$chr != "chrUn",]
intronsGR <- GRanges(seqnames = introns$chr,
                     ranges = IRanges(start = introns$start,
                                      end = introns$end),
                     strand = introns$strand)

# Get peaks and ranLoc that overlap CACTA_DTC
R1R3peaks_CACTA_DTC_overlap <- findOverlaps(query = CACTA_DTC_GR,
                                          subject = R1R3peaksGRselect,
                                          type = "any",
                                          select = "all",
                                          ignore.strand = TRUE)
CACTA_DTC_overlapping_R1R3peaksGR <- R1R3peaksGRselect[unique(subjectHits(R1R3peaks_CACTA_DTC_overlap))]
CACTA_DTC_overlapping_R1R3peaksGR_widths <- width(CACTA_DTC_overlapping_R1R3peaksGR)

R2Cpeaks_CACTA_DTC_overlap <- findOverlaps(query = CACTA_DTC_GR,
                                         subject = R2CpeaksGRselect,
                                         type = "any",
                                         select = "all",
                                         ignore.strand = TRUE)
CACTA_DTC_overlapping_R2CpeaksGR <- R2CpeaksGRselect[unique(subjectHits(R2Cpeaks_CACTA_DTC_overlap))]
CACTA_DTC_overlapping_R2CpeaksGR_widths <- width(CACTA_DTC_overlapping_R2CpeaksGR)

# Get peaks and ranLoc that overlap Copia_LTR_RLC
R1R3peaks_Copia_LTR_RLC_overlap <- findOverlaps(query = Copia_LTR_RLC_GR,
                                          subject = R1R3peaksGRselect,
                                          type = "any",
                                          select = "all",
                                          ignore.strand = TRUE)
Copia_LTR_RLC_overlapping_R1R3peaksGR <- R1R3peaksGRselect[unique(subjectHits(R1R3peaks_Copia_LTR_RLC_overlap))]
Copia_LTR_RLC_overlapping_R1R3peaksGR_widths <- width(Copia_LTR_RLC_overlapping_R1R3peaksGR)

R2Cpeaks_Copia_LTR_RLC_overlap <- findOverlaps(query = Copia_LTR_RLC_GR,
                                         subject = R2CpeaksGRselect,
                                         type = "any",
                                         select = "all",
                                         ignore.strand = TRUE)
Copia_LTR_RLC_overlapping_R2CpeaksGR <- R2CpeaksGRselect[unique(subjectHits(R2Cpeaks_Copia_LTR_RLC_overlap))]
Copia_LTR_RLC_overlapping_R2CpeaksGR_widths <- width(Copia_LTR_RLC_overlapping_R2CpeaksGR)

# Get peaks and ranLoc that overlap Gypsy_LTR_RLG
R1R3peaks_Gypsy_LTR_RLG_overlap <- findOverlaps(query = Gypsy_LTR_RLG_GR,
                                          subject = R1R3peaksGRselect,
                                          type = "any",
                                          select = "all",
                                          ignore.strand = TRUE)
Gypsy_LTR_RLG_overlapping_R1R3peaksGR <- R1R3peaksGRselect[unique(subjectHits(R1R3peaks_Gypsy_LTR_RLG_overlap))]
Gypsy_LTR_RLG_overlapping_R1R3peaksGR_widths <- width(Gypsy_LTR_RLG_overlapping_R1R3peaksGR)

R2Cpeaks_Gypsy_LTR_RLG_overlap <- findOverlaps(query = Gypsy_LTR_RLG_GR,
                                         subject = R2CpeaksGRselect,
                                         type = "any",
                                         select = "all",
                                         ignore.strand = TRUE)
Gypsy_LTR_RLG_overlapping_R2CpeaksGR <- R2CpeaksGRselect[unique(subjectHits(R2Cpeaks_Gypsy_LTR_RLG_overlap))]
Gypsy_LTR_RLG_overlapping_R2CpeaksGR_widths <- width(Gypsy_LTR_RLG_overlapping_R2CpeaksGR)




R1R3ranLoc_CACTA_DTC_overlap <- findOverlaps(query = CACTA_DTC_GR,
                                           subject = R1R3ranLocGR,
                                           type = "any",
                                           select = "all",
                                           ignore.strand = TRUE)
CACTA_DTC_overlapping_R1R3ranLocGR <- R1R3ranLocGRselect[unique(subjectHits(R1R3ranLoc_CACTA_DTC_overlap))]
CACTA_DTC_overlapping_R1R3ranLocGR_widths <- width(CACTA_DTC_overlapping_R1R3ranLocGR)

R2CranLoc_CACTA_DTC_overlap <- findOverlaps(query = CACTA_DTC_GR,
                                          subject = R2CranLocGR,
                                          type = "any",
                                          select = "all",
                                          ignore.strand = TRUE)
CACTA_DTC_overlapping_R2CranLocGR <- R2CranLocGR[unique(subjectHits(R2CranLoc_CACTA_DTC_overlap))]
CACTA_DTC_overlapping_R2CranLocGR_widths <- width(CACTA_DTC_overlapping_R2CranLocGR)


# Obtain 1000-bp gene promoters
# Get peaks that overlap gene promoters
promotersGR <- promoters(genesGR, upstream = 1000, downstream = 0)
R1R3peaks_promoters_overlap <- findOverlaps(query = promotersGR,
                                            subject = R1R3peaksGRselect,
                                            type = "any",
                                            select = "all",
                                            ignore.strand = TRUE)
promoter_overlapping_R1R3peaksGR <- R1R3peaksGR[unique(subjectHits(R1R3peaks_promoters_overlap))]
promoter_overlapping_R1R3peaksGR_widths <- width(promoter_overlapping_R1R3peaksGR)

R2Cpeaks_promoters_overlap <- findOverlaps(query = promotersGR,
                                           subject = R2CpeaksGRselect,
                                           type = "any",
                                           select = "all",
                                           ignore.strand = TRUE)
promoter_overlapping_R2CpeaksGR <- R2CpeaksGR[unique(subjectHits(R2Cpeaks_promoters_overlap))]
promoter_overlapping_R2CpeaksGR_widths <- width(promoter_overlapping_R2CpeaksGR)

R1R3ranLoc_promoters_overlap <- findOverlaps(query = promotersGR,
                                             subject = R1R3ranLocGR,
                                             type = "any",
                                             select = "all",
                                             ignore.strand = TRUE)
promoter_overlapping_R1R3ranLocGR <- R1R3ranLocGR[-unique(subjectHits(R1R3ranLoc_promoters_overlap))]
promoter_overlapping_R1R3ranLocGR_widths <- width(promoter_overlapping_R1R3ranLocGR)

R2CranLoc_promoters_overlap <- findOverlaps(query = promotersGR,
                                            subject = R2CranLocGR,
                                            type = "any",
                                            select = "all",
                                            ignore.strand = TRUE)
promoter_overlapping_R2CranLocGR <- R2CranLocGR[-unique(subjectHits(R2CranLoc_promoters_overlap))]
promoter_overlapping_R2CranLocGR_widths <- width(promoter_overlapping_R2CranLocGR)


# Get peaks that overlap gene regions immediately downstream of gene TSSs
# Obtain regions immediately downstream of gene TSSs (TSS to TSS+499 bp)
TSSsGR <- promoters(genesGR, upstream = 0, downstream = 500)

# Obtain regions immediately upstream of gene TTSs (TTS to TTS-499 bp)
source("/projects/ajt200/Rfunctions/TTSplus.R")
TTSsGR <- TTSplus(genesGR, upstream = 499, downstream = 0)
strand(TTSsGR) <- "*"
print(TTSsGR)

# Obtain 1000-bp gene terminators
terminatorsGR <- TTSplus(genesGR, upstream = -1, downstream = 1000)
strand(terminatorsGR) <- "*"
print(terminatorsGR)



# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = maskGR,
                                   subject = genesGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
genesGR <- genesGR[-subjectHits(mask_genes_overlap)]
# Retain strand information until after obtaining promoters, etc.

