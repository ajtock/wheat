#!/applications/R/R-3.5.0/bin/Rscript

# Run permutation tests evaluating wheat ASY1 peaks that overlap 
# H2A.Z peaks and H3K27me3 peaks for the relationship between
# recombination rate and overlap with loci matching a given motif
# that is enriched among those ASY1 peaks:

# Divide ASY1 peaks into those that do and do not overlap at least
# one motif-matching genomic locus.

# Define the subset of those peaks that do not overlap any motif-matching loci
# that have comparable base compositions to those that do.

# To control for regional differences in recombination rate, subset the peaks
# that do not overlap any motif-matching loci (and which have comparable base
# compositions to those that do) to include only those that are within 500 kb
# of at least one of those that do overlap a motif-matching locus.

# Compare mean 100-kb-scaled recombination rates (cM/Mb) for ASY1 peaks that
# do and do not overlap a motif-matching locus, using 10,000 random (permuted)
# sets of those that do not overlap a motif-matching locus to derive an
# empirical P-value.

# Usage:
#./motif_overlap_cMMb_permTest_region.R 100kb 1 ASY1_CS ASY1_CS_Rep1_ChIP 'p0.001_q0.01' 500000 10000 0.0001 'A' 'euchromatin' MAT1_TTCTTCTTCT 50

#winName <- "100kb"
#minMarkerDist <- "1"
#markChIP <- "ASY1_CS"
#libNameChIP <- "ASY1_CS_Rep1_ChIP"
#sigLevel <- "p0.001_q0.01"
#maxDistance <- 500000
#randomSets <- 10000
#minPval <- 0.0001
#genomeName <- "A"
#region <- "euchromatin"
#pwmName <- "MAT1_TTCTTCTTCT"
#size <- 50

args <- commandArgs(trailingOnly = T)
winName <- args[1]
minMarkerDist <- as.numeric(args[2])
markChIP <- args[3]
libNameChIP <- args[4]
sigLevel <- args[5]
maxDistance <- as.numeric(args[6])
randomSets <- as.numeric(args[7])
minPval <- as.numeric(args[8])
genomeName <- args[9]
region <- args[10]
pwmName <- args[11]
size <- as.numeric(args[12])

source("/projects/ajt200/Rfunctions/locus_midpoint.R")
library(Biostrings)
library(BSgenome.Taestivum.Cambridge.iwgscRefseqv1)
genome <- Taestivum
library(plotrix)

regionDir <- paste0(region, "/")
resDir <- paste0(regionDir, "results/")
plotDir <- paste0(regionDir, "histograms/")
gffDir <- "gff/"
bedDir <- "bed/"
system(paste0("[ -d ", regionDir, " ] || mkdir ", regionDir))
system(paste0("[ -d ", resDir, " ] || mkdir ", resDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
system(paste0("[ -d ", gffDir, " ] || mkdir ", gffDir))
system(paste0("[ -d ", bedDir, " ] || mkdir ", bedDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "pericentromeres") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}
# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "pericentromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Load position-weight matrix and obtain matches in the given subgenome partition
pwm <- as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/heterochromatin/motifs_ASY1_H2AZ_H3K27me3_peakOverlap_summits200bp/weeder2_bg_ranLoc_200bpseq/B/all/",
                                   pwmName, ".pwm"),
		    skip = 1, row.names = 1))
motif <- NULL
for(x in 1:dim(pwm)[2]) {
  mer <- names(pwm[,x][pwm[,x] == max(pwm[,x])])
  motif <- c(motif, mer)
}
motif <- paste0(motif, collapse = "")

# Extract genomeName-subgenome chromosomes
chrs <- chrs[grep(genomeName, chrs)]
genomeNums <- grep(genomeName,
                   names(genome))
chrList <- mclapply(genomeNums, function(i) {
  genome[[i]]
}, mc.cores = length(genomeNums))
names(chrList) <- names(genome)[genomeNums]

# Obtain motif matches in each genomeName-subgenome chromosome
hitsGR <- GRanges()
for(i in seq_along(chrList)) {
  print(chrs[i])
  hits <- matchPWM(pwm = pwm,
                   subject = chrList[[i]],
                   min.score = "99%")
  hitsGRchr <- GRanges(seqnames = chrs[i],
                       ranges = hits@ranges,
                       strand = "+")
  hitsRC <- matchPWM(pwm = reverseComplement(pwm),
                     subject = chrList[[i]],
                     min.score = "99%")
  hitsRCGRchr <- GRanges(seqnames = chrs[i],
                         ranges = hitsRC@ranges,
                         strand = "-")
  hitsGRchr <- sort(c(hitsGRchr, hitsRCGRchr))
  hitsGR <- append(hitsGR, hitsGRchr)
}
# Remove hits overlapping masked region
mask_hits_overlap <- findOverlaps(query = maskGR,
                                  subject = hitsGR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
hitsGR <- hitsGR[-subjectHits(mask_hits_overlap)]
print(hitsGR)
# Write as GRanges objects and in BED and GFF3 format
save(hitsGR,
     file = paste0(pwmName, "_motif_", motif, "_matchPWM_loci_in_",
                   genomeName, "genome_", region, ".RData"))
hitsgff <- data.frame(chr = as.character(seqnames(hitsGR)),
                      source = as.character(rep(".")),
                      feature = as.character(rep(paste0(pwmName[x], "_match"))),
                      start = as.integer(start(hitsGR)),
                      end = as.integer(end(hitsGR)),
                      score = as.character(rep(".")),
                      strand = as.character(strand(hitsGR)),
                      phase = as.character(rep(".")),
                      attribute = as.character(rep(".")))
write.table(hitsgff,
            file = paste0(gffDir, pwmName, "_motif_", motif, "_matchPWM_loci_in_",
                          genomeName, "genome_", region, ".gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")
hitsbed <- data.frame(chr = as.character(seqnames(hitsGR)),
                      start = as.integer(start(hitsGR)-1),
                      end = as.integer(end(hitsGR)),
                      name = as.integer(1:length(hitsGR)),
                      score = rep("NA", length(hitsGR)),
                      strand = as.character(strand(hitsGR)))
write.table(hitsbed,
            file = paste0(bedDir, pwmName, "_motif_", motif, "_matchPWM_loci_in_",
                          genomeName, "genome_", region, ".bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")


# Load ASY1 peaks that overlap H2A.Z and H3K27me3 peaks,
# and re-define as summits flanked by size/2 bp
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/", markChIP, "/snakemake_ChIPseq/mapped/",
                           "ASY1peakProfiles_overlap_H2AZ_H3K27me3_peaks/heatmaps/gff/bodies_by_log2_",
                           libNameChIP, "_control/",
                           libNameChIP, "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel, "_noMinWidth_in_",
                           genomeName, "genome_", region,
                           "_overlapping_H2AZ_Rep1_ChIP_and_H3K27me3_ChIP_SRR6350666_peaks",
                           "_ordered_by_log2_", libNameChIP, "_control_in_bodies.gff"))
colnames(peaks) <- c("chr", "source", "feature",
                     "start", "end", "sigval",
                     "strand", "qval", "summit0based")
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start,
                                    end = peaks$end),
                   strand = "*",
                   sigval = peaks$sigval,
                   qval = peaks$qval,
		   summit0based = peaks$summit0based)
peaksGR_summits <- GRanges(seqnames = seqnames(peaksGR),
                           ranges = IRanges(start = start(peaksGR)+peaksGR$summit0based,
                                            end = start(peaksGR)+peaksGR$summit0based),
                           strand = "*",
                           sigval = peaksGR$sigval,
                           qval = peaksGR$qval,
                           summit0based = peaksGR$summit0based)
peaksGR_summitsFlank <- locMidpointFlank(x = peaksGR_summits,
                                         leftFlank = (size/2),
                                         rightFlank = (size/2))

# Convert windowed recombination rate into GRanges
cMMb <- read.table(paste0(
                   "/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/",
                   "cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                   as.character(minMarkerDist), "bp_", winName, ".txt"))
cMMbGR <- GRanges(seqnames = cMMb$chr,
                  ranges = IRanges(start = cMMb$windowStart,
                                   end = cMMb$windowEnd),
                  strand = "*",
                  cMMb = cMMb$cMMb)

# Obtain winName-scaled cMMb values for each peak
# Where features overlap more than one winName window, calculate mean cMMb
peaks_cMMb_overlap <- findOverlaps(query = peaksGR_summitsFlank,
                                   subject = cMMbGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaks_cMMb_overlapList <- lapply(seq_along(peaksGR_summitsFlank), function(x) {
  subjectHits(peaks_cMMb_overlap)[queryHits(peaks_cMMb_overlap) == x]
})
peaks_cMMb <- sapply(peaks_cMMb_overlapList,
                     function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
peaksGR <- GRanges(peaksGR_summitsFlank,
                   sigval = peaksGR_summitsFlank$sigval,
                   qval = peaksGR_summitsFlank$qval,
                   summit0based = peaksGR_summitsFlank$summit0based,
                   cMMb = peaks_cMMb)

# Obtain peaks that do and do not overlap motif-matching loci
hits_peaks_overlap <- findOverlaps(query = hitsGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaksGR_hitsGR <- peaksGR[unique(subjectHits(hits_peaks_overlap))]
peaksGR_no_hitsGR <- peaksGR[-subjectHits(hits_peaks_overlap)]
peaksGR_hitsGR_extend <- GRanges(seqnames = seqnames(peaksGR_hitsGR),
                                 ranges = IRanges(start = start(peaksGR_hitsGR)-maxDistance,
                                                  end = end(peaksGR_hitsGR)+maxDistance),
                                 strand = strand(peaksGR_hitsGR),
                                 cMMb = peaksGR_hitsGR$cMMb)
peaksGR_no_hitsGR_nearby <- findOverlaps(query = peaksGR_hitsGR_extend,
                                         subject = peaksGR_no_hitsGR,
                                         type = "any",
                                         select = "all",
                                         ignore.strand = TRUE)

peaksGR_no_hitsGR <- peaksGR_no_hitsGR[unique(subjectHits(peaksGR_no_hitsGR_nearby))]

print(peaksGR_hitsGR)
print(peaksGR_no_hitsGR)
stopifnot(length(peaksGR_hitsGR) +
          length(peaksGR_no_hitsGR) == dim(peaks)[1])
# Obtain sequences
peaksGR_hitsGRseq <- DNAStringSet()
peaksGR_no_hitsGRseq <- DNAStringSet()
for(h in seq_along(chrs)) {
  peaksGR_hitsGRseq_chr <- DNAStringSet()
  peaksGR_no_hitsGRseq_chr <- DNAStringSet()
  # peaks overlapping motif-matching loci
  for(i in seq_along(peaksGR_hitsGR[seqnames(peaksGR_hitsGR) == chrs[h]])) {
    peaksGR_hitsGRseq_chr <- c(peaksGR_hitsGRseq_chr,
                               DNAStringSet(chrList[[h]][(start(peaksGR_hitsGR[seqnames(peaksGR_hitsGR)
                                                                               == chrs[h]][i])):
                                                         (end(peaksGR_hitsGR[seqnames(peaksGR_hitsGR)
                                                                             == chrs[h]][i]))]))
  }
  peaksGR_hitsGRseq <- c(peaksGR_hitsGRseq, peaksGR_hitsGRseq_chr)
  # peaks not overlapping motif-matching loci
  for(i in seq_along(peaksGR_no_hitsGR[seqnames(peaksGR_no_hitsGR) == chrs[h]])) {
    peaksGR_no_hitsGRseq_chr <- c(peaksGR_no_hitsGRseq_chr,
                                  DNAStringSet(chrList[[h]][(start(peaksGR_no_hitsGR[seqnames(peaksGR_no_hitsGR)
                                                                                     == chrs[h]][i])):
                                                            (end(peaksGR_no_hitsGR[seqnames(peaksGR_no_hitsGR)
                                                                                   == chrs[h]][i]))]))
  }
  peaksGR_no_hitsGRseq <- c(peaksGR_no_hitsGRseq, peaksGR_no_hitsGRseq_chr)
}

# Obtain base frequencies over central 100-bp sequences
# peaks overlapping motif-matching loci
peaksGR_hitsGRseq_AtoN <- NULL
for(i in 1:length(peaksGR_hitsGRseq)) {
  peaksGR_hitsGRseq_AtoN <- rbind(peaksGR_hitsGRseq_AtoN, 
                                  alphabetFrequency(peaksGR_hitsGRseq[[i]][
                                                      ((size/2)-49):((size/2)+51)
                                                    ])[c(1:4,15)])
}
peaksGR_hitsGRseq_AtoN <- as.data.frame(peaksGR_hitsGRseq_AtoN)
peaksGR_hitsGRseq_AtoN_ranges <- data.frame(min = c(min(peaksGR_hitsGRseq_AtoN$A),
                                                    min(peaksGR_hitsGRseq_AtoN$C),
                                                    min(peaksGR_hitsGRseq_AtoN$G),
                                                    min(peaksGR_hitsGRseq_AtoN$T),
                                                    min(peaksGR_hitsGRseq_AtoN$N)),
                                            max = c(max(peaksGR_hitsGRseq_AtoN$A),
                                                    max(peaksGR_hitsGRseq_AtoN$C),
                                                    max(peaksGR_hitsGRseq_AtoN$G),
                                                    max(peaksGR_hitsGRseq_AtoN$T),
                                                    max(peaksGR_hitsGRseq_AtoN$N))) 
rownames(peaksGR_hitsGRseq_AtoN_ranges) <- c("A", "C", "G", "T", "N")
# peaks not overlapping motif-matching loci
peaksGR_no_hitsGRseq_AtoN <- NULL
for(i in 1:length(peaksGR_no_hitsGRseq)) {
  peaksGR_no_hitsGRseq_AtoN <- rbind(peaksGR_no_hitsGRseq_AtoN,
                                     alphabetFrequency(peaksGR_no_hitsGRseq[[i]][
                                                         ((size/2)-49):((size/2)+51)
                                                       ])[c(1:4,15)])
}
peaksGR_no_hitsGRseq_AtoN <- as.data.frame(peaksGR_no_hitsGRseq_AtoN)
# select those peaks that do not overlap motif-matching loci
# with base frequencies within the ranges of those that do
peaksGR_no_hitsGRseq
test <- peaksGR_no_hitsGRseq_AtoN[# A
                                  peaksGR_no_hitsGRseq_AtoN$A >=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "A",]$min &
                                  peaksGR_no_hitsGRseq_AtoN$A <=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "A",]$max &
                                  # C
                                  peaksGR_no_hitsGRseq_AtoN$C >=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "C",]$min &
                                  peaksGR_no_hitsGRseq_AtoN$C <=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "C",]$max &
                                  # G
                                  peaksGR_no_hitsGRseq_AtoN$G >=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "G",]$min &
                                  peaksGR_no_hitsGRseq_AtoN$G <=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "G",]$max &
                                  # T
                                  peaksGR_no_hitsGRseq_AtoN$T >=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "T",]$min &
                                  peaksGR_no_hitsGRseq_AtoN$T <=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "T",]$max &
                                  # N
                                  peaksGR_no_hitsGRseq_AtoN$N >=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "N",]$min &
                                  peaksGR_no_hitsGRseq_AtoN$N <=
                                  peaksGR_hitsGRseq_AtoN_ranges[
                                    rownames(peaksGR_hitsGRseq_AtoN_ranges) == "N",]$max,
]

peaksGR_no_hitsGRseq_AtoN_ranges <- data.frame(min = c(min(peaksGR_no_hitsGRseq_AtoN$A),
                                                       min(peaksGR_no_hitsGRseq_AtoN$C),
                                                       min(peaksGR_no_hitsGRseq_AtoN$G),
                                                       min(peaksGR_no_hitsGRseq_AtoN$T),
                                                       min(peaksGR_no_hitsGRseq_AtoN$N)),
                                               max = c(max(peaksGR_no_hitsGRseq_AtoN$A),
                                                       max(peaksGR_no_hitsGRseq_AtoN$C),
                                                       max(peaksGR_no_hitsGRseq_AtoN$G),
                                                       max(peaksGR_no_hitsGRseq_AtoN$T),
                                                       max(peaksGR_no_hitsGRseq_AtoN$N)))
rownames(peaksGR_no_hitsGRseq_AtoN_ranges) <- c("A", "C", "G", "T", "N")




# Load table of representative genes and convert into GRanges
genes <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA.gff3",
                    colClasses = c(NA,
                                   rep("NULL", 2),
                                   rep(NA, 2),
                                   "NULL", NA, "NULL", NA))
colnames(genes) <- c("chr", "start", "end", "strand", "geneID")
genes <- genes[genes$chr != "chrUn",]
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = genes$strand,
                   geneID = genes$geneID)
print(genesGR)
# Subset to include only those overlapping specified region (e.g., euchromatin)
genesGR <- subsetByOverlaps(x = genesGR,
                            ranges = regionGR,
                            type = "any",
                            invert = FALSE,
                            ignore.strand = TRUE)

# Obtain 1000-bp gene promoters and terminators
promotersGR <- promoters(genesGR, upstream = 1000, downstream = 0)
print(promotersGR)

# Obtain regions immediately downstream of gene TSSs (TSS to TSS+499 bp)
TSSsGR <- promoters(genesGR, upstream = 0, downstream = 500)
print(TSSsGR)

# Obtain regions immediately upstream of gene TTSs (TTS to TTS-499 bp)
source("/projects/ajt200/Rfunctions/TTSplus.R")
TTSsGR <- TTSplus(genesGR, upstream = 499, downstream = 0)
print(TTSsGR)

terminatorsGR <- TTSplus(genesGR, upstream = -1, downstream = 1000)
print(terminatorsGR)


# Convert windowed recombination rate into GRanges
cMMb <- read.table(paste0(
                   "/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/",
                   "cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                   as.character(minMarkerDist), "bp_", winName, ".txt"))
cMMbGR <- GRanges(seqnames = cMMb$chr,
                  ranges = IRanges(start = cMMb$windowStart,
                                   end = cMMb$windowEnd),
                  strand = "*",
                  cMMb = cMMb$cMMb)

# Obtain winName-scaled cMMb values for each gene promoter and terminator
# Where features overlap more than one winName window, calculate mean cMMb
# genes
gene_cMMb_overlaps <- findOverlaps(query = genesGR,
                                   subject = cMMbGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
gene_cMMb_overlapsList <- lapply(seq_along(genesGR), function(x) {
  subjectHits(gene_cMMb_overlaps)[queryHits(gene_cMMb_overlaps) == x]
})
gene_cMMb <- sapply(gene_cMMb_overlapsList,
                    function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
genesGR <- GRanges(genesGR,
                   geneID = genesGR$geneID,
                   cMMb = gene_cMMb)

# promoters
promoter_cMMb_overlaps <- findOverlaps(query = promotersGR,
                                       subject = cMMbGR,
                                       type = "any",
                                       select = "all",
                                       ignore.strand = TRUE)
promoter_cMMb_overlapsList <- lapply(seq_along(promotersGR), function(x) {
  subjectHits(promoter_cMMb_overlaps)[queryHits(promoter_cMMb_overlaps) == x]
})
## OR
#promoter_cMMb_overlapsList <- getOverlaps(coordinates = promotersGR,
#                                          segments = cMMbGR,
#                                          overlapType = "overlapping",
#                                          whichOverlaps = TRUE,
#                                          ignoreStrand = TRUE)
promoter_cMMb <- sapply(promoter_cMMb_overlapsList,
                        function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
promotersGR <- GRanges(promotersGR,
                       geneID = promotersGR$geneID,
                       cMMb = promoter_cMMb)

# TSSs
TSS_cMMb_overlaps <- findOverlaps(query = TSSsGR,
                                  subject = cMMbGR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
TSS_cMMb_overlapsList <- lapply(seq_along(TSSsGR), function(x) {
  subjectHits(TSS_cMMb_overlaps)[queryHits(TSS_cMMb_overlaps) == x]
})
TSS_cMMb <- sapply(TSS_cMMb_overlapsList,
                   function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
TSSsGR <- GRanges(TSSsGR,
                  geneID = TSSsGR$geneID,
                  cMMb = TSS_cMMb)

# TTSs
TTS_cMMb_overlaps <- findOverlaps(query = TTSsGR,
                                  subject = cMMbGR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
TTS_cMMb_overlapsList <- lapply(seq_along(TTSsGR), function(x) {
  subjectHits(TTS_cMMb_overlaps)[queryHits(TTS_cMMb_overlaps) == x]
})
TTS_cMMb <- sapply(TTS_cMMb_overlapsList,
                   function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
TTSsGR <- GRanges(TTSsGR,
                  geneID = TTSsGR$geneID,
                  cMMb = TTS_cMMb)

# terminators
terminator_cMMb_overlaps <- findOverlaps(query = terminatorsGR,
                                         subject = cMMbGR,
                                         type = "any",
                                         select = "all",
                                         ignore.strand = TRUE)
terminator_cMMb_overlapsList <- lapply(seq_along(terminatorsGR), function(x) {
  subjectHits(terminator_cMMb_overlaps)[queryHits(terminator_cMMb_overlaps) == x]
})
terminator_cMMb <- sapply(terminator_cMMb_overlapsList,
                          function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
terminatorsGR <- GRanges(terminatorsGR,
                         geneID = terminatorsGR$geneID,
                         cMMb = terminator_cMMb)

# Load peaks for given ChIP-seq dataset
if(libNameChIP %in% c("H3K4me3_ChIP_SRR6350668",
                      "H3K27me3_ChIP_SRR6350666",
                      "H3K36me3_ChIP_SRR6350670",
                      "H3K9ac_ChIP_SRR6350667",
                      "CENH3_ChIP_SRR1686799")) {
  peaksFile <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                      markChIP, "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                      sigLevel, "/",
                      libNameChIP, "_rangerPeaksGRmergedOverlaps_minuslog10_",
                      sigLevel, "_noMinWidth.RData")
} else {
  peaksFile <- paste0("/home/ajt200/analysis/wheat/",
                      markChIP, "/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/",
                      sigLevel, "/",
                      libNameChIP, "_rangerPeaksGRmergedOverlaps_minuslog10_",
                      sigLevel, "_noMinWidth.RData")
}

load(peaksFile)
peaksGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
# Subset to include only those overlapping specified region (e.g., euchromatin)
peaksGR <- subsetByOverlaps(x = peaksGR,
                            ranges = regionGR,
                            type = "any",
                            invert = FALSE,
                            ignore.strand = TRUE)

# Obtain gene parts that do or don't overlap peaks
# genes
gene_peak_overlapsGR <- subsetByOverlaps(x = genesGR,
                                         ranges = peaksGR,
                                         type = "any",
                                         invert = FALSE,
                                         ignore.strand = TRUE)
gene_peak_nonoverlapsGR <- subsetByOverlaps(x = genesGR,
                                            ranges = peaksGR,
                                            type = "any",
                                            invert = TRUE,
                                            ignore.strand = TRUE)
gene_peak_overlapsGR_extend <- GRanges(seqnames = seqnames(gene_peak_overlapsGR),
                                       ranges = IRanges(start = start(gene_peak_overlapsGR)-maxDistance,
                                                        end = end(gene_peak_overlapsGR)+maxDistance),
                                       strand = strand(gene_peak_overlapsGR),
                                       cMMb = gene_peak_overlapsGR$cMMb)
gene_peak_nonoverlapsGR_nearby <- subsetByOverlaps(x = gene_peak_nonoverlapsGR,
                                                   ranges = gene_peak_overlapsGR_extend,
                                                   type = "any",
                                                   invert = FALSE,
                                                   ignore.strand = TRUE)
# promoters
promoter_peak_overlapsGR <- subsetByOverlaps(x = promotersGR,
                                             ranges = peaksGR,
                                             type = "any",
                                             invert = FALSE,
                                             ignore.strand = TRUE)
promoter_peak_nonoverlapsGR <- subsetByOverlaps(x = promotersGR,
                                                ranges = peaksGR,
                                                type = "any",
                                                invert = TRUE,
                                                ignore.strand = TRUE)
promoter_peak_overlapsGR_extend <- GRanges(seqnames = seqnames(promoter_peak_overlapsGR),
                                           ranges = IRanges(start = start(promoter_peak_overlapsGR)-maxDistance,
                                                            end = end(promoter_peak_overlapsGR)+maxDistance),
                                           strand = strand(promoter_peak_overlapsGR),
                                           cMMb = promoter_peak_overlapsGR$cMMb)
promoter_peak_nonoverlapsGR_nearby <- subsetByOverlaps(x = promoter_peak_nonoverlapsGR,
                                                       ranges = promoter_peak_overlapsGR_extend,
                                                       type = "any",
                                                       invert = FALSE,
                                                       ignore.strand = TRUE)
# TSSs
TSS_peak_overlapsGR <- subsetByOverlaps(x = TSSsGR,
                                        ranges = peaksGR,
                                        type = "any",
                                        invert = FALSE,
                                        ignore.strand = TRUE)
TSS_peak_nonoverlapsGR <- subsetByOverlaps(x = TSSsGR,
                                           ranges = peaksGR,
                                           type = "any",
                                           invert = TRUE,
                                           ignore.strand = TRUE)
TSS_peak_overlapsGR_extend <- GRanges(seqnames = seqnames(TSS_peak_overlapsGR),
                                      ranges = IRanges(start = start(TSS_peak_overlapsGR)-maxDistance,
                                                       end = end(TSS_peak_overlapsGR)+maxDistance),
                                      strand = strand(TSS_peak_overlapsGR),
                                      cMMb = TSS_peak_overlapsGR$cMMb)
TSS_peak_nonoverlapsGR_nearby <- subsetByOverlaps(x = TSS_peak_nonoverlapsGR,
                                                  ranges = TSS_peak_overlapsGR_extend,
                                                  type = "any",
                                                  invert = FALSE,
                                                  ignore.strand = TRUE)
# TTSs
TTS_peak_overlapsGR <- subsetByOverlaps(x = TTSsGR,
                                        ranges = peaksGR,
                                        type = "any",
                                        invert = FALSE,
                                        ignore.strand = TRUE)
TTS_peak_nonoverlapsGR <- subsetByOverlaps(x = TTSsGR,
                                           ranges = peaksGR,
                                           type = "any",
                                           invert = TRUE,
                                           ignore.strand = TRUE)
TTS_peak_overlapsGR_extend <- GRanges(seqnames = seqnames(TTS_peak_overlapsGR),
                                      ranges = IRanges(start = start(TTS_peak_overlapsGR)-maxDistance,
                                                       end = end(TTS_peak_overlapsGR)+maxDistance),
                                      strand = strand(TTS_peak_overlapsGR),
                                      cMMb = TTS_peak_overlapsGR$cMMb)
TTS_peak_nonoverlapsGR_nearby <- subsetByOverlaps(x = TTS_peak_nonoverlapsGR,
                                                  ranges = TTS_peak_overlapsGR_extend,
                                                  type = "any",
                                                  invert = FALSE,
                                                  ignore.strand = TRUE)
# terminators
terminator_peak_overlapsGR <- subsetByOverlaps(x = terminatorsGR,
                                               ranges = peaksGR,
                                               type = "any",
                                               invert = FALSE,
                                               ignore.strand = TRUE)
terminator_peak_nonoverlapsGR <- subsetByOverlaps(x = terminatorsGR,
                                                  ranges = peaksGR,
                                                  type = "any",
                                                  invert = TRUE,
                                                  ignore.strand = TRUE)
terminator_peak_overlapsGR_extend <- GRanges(seqnames = seqnames(terminator_peak_overlapsGR),
                                             ranges = IRanges(start = start(terminator_peak_overlapsGR)-maxDistance,
                                                              end = end(terminator_peak_overlapsGR)+maxDistance),
                                             strand = strand(terminator_peak_overlapsGR),
                                             cMMb = terminator_peak_overlapsGR$cMMb)
terminator_peak_nonoverlapsGR_nearby <- subsetByOverlaps(x = terminator_peak_nonoverlapsGR,
                                                         ranges = terminator_peak_overlapsGR_extend,
                                                         type = "any",
                                                         invert = FALSE,
                                                         ignore.strand = TRUE)

# Define function to randomly select n rows from
# a GRanges object (features)
selectRandomFeatures <- function(features, n) {
  return(features[sample(x = length(features),
                         size = n,
                         replace = FALSE)])
}

# Define seed so that random selections are reproducible
set.seed(93750174)

# Set class for permutation test results object
setClass("permTest_cMMb",
         representation(alternative = "character",
                        alpha0.05 = "numeric",
                        pval = "numeric",
                        observed = "numeric",
                        permuted = "numeric",
                        expected = "numeric"))

# Disable scientific notation (e.g., 0.0001 rather than 1e-04)
options(scipen = 100)

# Permutation test function to evaluate if mean cM/Mb at peak-containing loci is
# significantly higher or lower than at random sets of peak-less loci
loci_cMMb_permTest <- function(targets,
                               nontargets,
                               targetsName,
                               targetsNamePlot,
                               genomeName,
                               genomeNamePlot,
                               resultsDir,
                               plotDir) {
  # Select genome-specific chromosomes (e.g., in the "A" genome)
  chrs <- chrs[grep(genomeName, chrs)]
   
  # Apply selectRandomFeatures() function on a per-chromosome basis
  # and append the selected ranges to a growing GRanges object (ranLocGR)
  # Repeat randomSets times to create a GRangesList object (ranLocGRL)
  ranLocGRL <- mclapply(1:randomSets, function(x) {
    ranLocGR <- GRanges()
    for(i in 1:length(chrs)) {
      ranLocChrGR <- selectRandomFeatures(features = nontargets[seqnames(nontargets)
                                                                == chrs[i]],
                                          n = length(targets[seqnames(targets)
                                                             == chrs[i]]))
      ranLocGR <- append(ranLocGR, ranLocChrGR)
    }
    ranLocGR
  }, mc.cores = 6)
  
  # Calculate mean cM/Mb for peak-containing loci (targets) and for
  # each set of peak-less random loci (selected from nontargets)
  targets_cMMbMean <- mean(targets$cMMb, na.rm = T)
  ranLoc_cMMbMean <- unlist(mclapply(seq_along(ranLocGRL), function(x) {
    mean(ranLocGRL[[x]]$cMMb, na.rm = T)
  }, mc.cores = 6))
  
  # Determine whether mean cM/Mb values at peak-less random loci are lower than or
  # higher than at peak-containing loci
  ranLoc_cMMb_lessThan_targets_cMMb_Bool <- sapply(seq_along(ranLoc_cMMbMean),
    function(x) {
      ranLoc_cMMbMean[x] < targets_cMMbMean
  })
  ranLoc_cMMb_moreThan_targets_cMMb_Bool <- sapply(seq_along(ranLoc_cMMbMean),
    function(x) {
      ranLoc_cMMbMean[x] > targets_cMMbMean
  })
  
  # Calculate P-values and significance levels
  if(mean(ranLoc_cMMbMean) < targets_cMMbMean) {
  #if(sum(ranLoc_cMMb_lessThan_targets_cMMb_Bool)
  #   > (length(ranLoc_cMMb_lessThan_targets_cMMb_Bool)/2)) {
    pval <- 1-(sum(ranLoc_cMMb_lessThan_targets_cMMb_Bool)/randomSets)
    if(pval == 0) {
      pval <- minPval
    }
    MoreOrLessThanRandom <- "MoreThanRandom"
    alpha0.05 <- quantile(ranLoc_cMMbMean, 0.95)[[1]]
  } else {
    pval <- 1-(sum(ranLoc_cMMb_moreThan_targets_cMMb_Bool)/randomSets)
    if(pval == 0) {
      pval <- minPval
    }
    MoreOrLessThanRandom <- "LessThanRandom"
    alpha0.05 <- quantile(ranLoc_cMMbMean, 0.05)[[1]]
  }
   
  # Create permtation test results object
  permTestResults <- new("permTest_cMMb",
                         alternative = MoreOrLessThanRandom,
                         alpha0.05 = alpha0.05,
                         pval = pval,
                         observed = targets_cMMbMean,
                         permuted = ranLoc_cMMbMean,
                         expected = mean(ranLoc_cMMbMean))
  save(permTestResults,
       file = paste0(resultsDir,
                     targetsName, "_vs_",
                     as.character(randomSets), "randomSets_",
                     "minInterMarkerDist", as.character(minMarkerDist), "bp_",
                     winName, "Scaled_cMMb_permTestResults.RData"))
  
  # Generate histogram
  pdf(paste0(plotDir,
             "hist_",
             targetsName, "_vs_",
             as.character(randomSets), "randomSets_",
             "minInterMarkerDist", as.character(minMarkerDist), "bp_",
             winName, "Scaled_cMMb_permTestResults.pdf"),
             height = 4.5, width = 5)
  par(mar = c(3.1, 3.1, 4.1, 1.1),
      mgp = c(1.85, 0.75, 0))
  ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
  #options(scipen = 100)
  # Calculate max density
  maxDensityPlus <- max(density(permTestResults@permuted)$y)*1.2
  if(permTestResults@alternative == "MoreThanRandom") {
    xlim <- c(pmax(0, min(permTestResults@permuted)/1.2),
              pmax(permTestResults@observed*1.2, permTestResults@alpha0.05*1.2))
    textX1 <- pmax(0.05, min(permTestResults@permuted)/1.15)
  } else {
    xlim <- c(pmax(0, permTestResults@observed/1.2),
              max(permTestResults@permuted)*1.2)
    textX1 <- max(permTestResults@permuted)*1.15
  }
  hist(permTestResults@permuted,
       freq = FALSE,
       col = "dodgerblue",
       border = NA,
       lwd = 2,
       xlim = c(pretty(xlim)[1],
                pretty(xlim)[length(pretty(xlim))]),
       ylim = c(0,
                maxDensityPlus),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "",
       axes = FALSE)
  axis(side = 2,
       at = pretty(density(permTestResults@permuted)$y),
       lwd = 2)
  mtext(side = 2,
        text = "Density",
        line = 1.85)
  axis(side = 1,
       at = pretty(xlim),
       lwd = 2)
  mtext(side = 1,
        text = "Mean cM/Mb",
        line = 1.85)
  titleText <- list(bquote(.(as.character(targetsNamePlot))),
                    bquote(italic("P")*" = "*
                           .(as.character(round(permTestResults@pval,
                                                digits = 6)))),
                    bquote("Permutations = "*.(prettyNum(randomSets,
                                                         big.mark = ",",
                                                         trim = T))))
  mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = c(0.7, 1, 1))
  lines(density(permTestResults@permuted),
        col = "dodgerblue3",
        lwd = 1.5)
  ablineclip(v = permTestResults@expected,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
  ablineclip(v = permTestResults@observed,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
  ablineclip(v = permTestResults@alpha0.05,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
  text(x = c(textX1,
             permTestResults@expected,
             permTestResults@observed,
             permTestResults@alpha0.05),
       y = c(maxDensityPlus*.95,
             maxDensityPlus,
             maxDensityPlus,
             maxDensityPlus*.95),
       labels = c("Permuted",
                  "Expected",
                  "Observed",
                  expression(alpha*" = 0.05")),
       col = c("dodgerblue",
               "black",
               "forestgreen",
               "red"),
       cex = 0.8)
  dev.off()

}

# Apply to promoters and terminators in each genome
genomeNames <- c("A", "B", "D", "")
genomeNamesPlot <- c("A", "B", "D", "all")
sapply(seq_along(genomeNames), function(z) {
  loci_cMMb_permTest(targets = gene_peak_overlapsGR[grep(genomeNames[z],
                                                         seqnames(gene_peak_overlapsGR))],
                     nontargets = gene_peak_nonoverlapsGR_nearby[grep(genomeNames[z],
                                                                      seqnames(gene_peak_nonoverlapsGR_nearby))],
                     targetsName = paste0("genes_in_", genomeNames[z], "genome_", region,
                                          "_overlapping_", libNameChIP, "_peaks"),
                     targetsNamePlot = paste0("Genes in ", genomeNamesPlot[z], "-genome ", region,
                                              " overlapping ", libNameChIP, " peaks"),
                     genomeName = genomeNames[z],
                     genomeNamePlot = genomeNamesPlot[z],
                     resultsDir = geneResDir,
                     plotDir = genePlotDir)
  loci_cMMb_permTest(targets = promoter_peak_overlapsGR[grep(genomeNames[z],
                                                             seqnames(promoter_peak_overlapsGR))],
                     nontargets = promoter_peak_nonoverlapsGR_nearby[grep(genomeNames[z],
                                                                          seqnames(promoter_peak_nonoverlapsGR_nearby))],
                     targetsName = paste0("gene_promoters_in_", genomeNames[z], "genome_", region,
                                          "_overlapping_", libNameChIP, "_peaks"),
                     targetsNamePlot = paste0("Gene promoters in ", genomeNamesPlot[z], "-genome ", region,
                                              " overlapping ", libNameChIP, " peaks"),
                     genomeName = genomeNames[z],
                     genomeNamePlot = genomeNamesPlot[z],
                     resultsDir = promoterResDir,
                     plotDir = promoterPlotDir)
  loci_cMMb_permTest(targets = TSS_peak_overlapsGR[grep(genomeNames[z],
                                                        seqnames(TSS_peak_overlapsGR))],
                     nontargets = TSS_peak_nonoverlapsGR_nearby[grep(genomeNames[z],
                                                                     seqnames(TSS_peak_nonoverlapsGR_nearby))],
                     targetsName = paste0("gene_TSSs_in_", genomeNames[z], "genome_", region,
                                          "_overlapping_", libNameChIP, "_peaks"),
                     targetsNamePlot = paste0("Gene TSSs+500 bp in ", genomeNamesPlot[z], "-genome ", region,
                                              " overlapping ", libNameChIP, " peaks"),
                     genomeName = genomeNames[z],
                     genomeNamePlot = genomeNamesPlot[z],
                     resultsDir = TSSResDir,
                     plotDir = TSSPlotDir)
  loci_cMMb_permTest(targets = TTS_peak_overlapsGR[grep(genomeNames[z],
                                                        seqnames(TTS_peak_overlapsGR))],
                     nontargets = TTS_peak_nonoverlapsGR_nearby[grep(genomeNames[z],
                                                                     seqnames(TTS_peak_nonoverlapsGR_nearby))],
                     targetsName = paste0("gene_TTSs_in_", genomeNames[z], "genome_", region,
                                          "_overlapping_", libNameChIP, "_peaks"),
                     targetsNamePlot = paste0("Gene TTSs-500 bp in ", genomeNamesPlot[z], "-genome ", region,
                                              " overlapping ", libNameChIP, " peaks"),
                     genomeName = genomeNames[z],
                     genomeNamePlot = genomeNamesPlot[z],
                     resultsDir = TTSResDir,
                     plotDir = TTSPlotDir)
  loci_cMMb_permTest(targets = terminator_peak_overlapsGR[grep(genomeNames[z],
                                                               seqnames(terminator_peak_overlapsGR))],
                     nontargets = terminator_peak_nonoverlapsGR_nearby[grep(genomeNames[z],
                                                                            seqnames(terminator_peak_nonoverlapsGR_nearby))],
                     targetsName = paste0("gene_terminators_in_", genomeNames[z], "genome_", region,
                                          "_overlapping_", libNameChIP, "_peaks"),
                     targetsNamePlot = paste0("Gene terminators in ", genomeNamesPlot[z], "-genome ", region,
                                              " overlapping ", libNameChIP, " peaks"),
                     genomeName = genomeNames[z],
                     genomeNamePlot = genomeNamesPlot[z],
                     resultsDir = terminatorResDir,
                     plotDir = terminatorPlotDir)
})
