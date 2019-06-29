#!/applications/R/R-3.5.0/bin/Rscript

# Run permutation tests evaluating wheat gene promoters (1-kb region upstream of TSS),
# regions immediately downstream of gene TSSs (TSS to TSS+499 bp),
# regions immediately upstream of gene TTSs (TTS to TTS-499 bp), and
# gene terminators (1-kb downstream of TTS) for the relationship between
# recombination rate and overlap with peaks in a given ChIP-seq dataset:

# Divide genes into those whose parts (e.g., promoters) overlap or do not overlap
# at least one peak.

# To control for regional differences in recombination rate, subset gene parts that
# do not overlap a peak to include only those that are within 500 kb of at least one
# of those that do overlap a peak.

# Compare mean 100-kb-scaled recombination rates (cM/Mb) for gene parts that do and
# do not overlap a peak, using 10,000 random (permuted) sets of those that do not
# overlap a peak to derive an empirical P-value.

# Usage:
#./peak_overlap_cMMb_permTest_region.R 100kb 1 H3K4me3 H3K4me3_Rep1_ChIP 'p0.05_q0.05' 500000 10000 0.0001 'euchromatin'

#winName <- "100kb"
#minMarkerDist <- "1"
#markChIP <- "H3K4me3"
#libNameChIP <- "H3K4me3_Rep1_ChIP"
#sigLevel <- "p0.05_q0.05"
#maxDistance <- 500000
#randomSets <- 10000
#minPval <- 0.0001
#region <- "euchromatin"

args <- commandArgs(trailingOnly = T)
winName <- args[1]
minMarkerDist <- as.numeric(args[2])
markChIP <- args[3]
libNameChIP <- args[4]
sigLevel <- args[5]
maxDistance <- as.numeric(args[6])
randomSets <- as.numeric(args[7])
minPval <- as.numeric(args[8])
region <- args[9]

library(GenomicRanges)
library(parallel)
library(plotrix)

regionDir <- paste0(region, "/")
geneDir <- paste0(regionDir, "genes/")
promoterDir <- paste0(regionDir, "promoters/")
TSSDir <- paste0(regionDir, "TSSplus500bp/")
TTSDir <- paste0(regionDir, "TTSminus500bp/")
terminatorDir <- paste0(regionDir, "terminators/")
geneResDir <- paste0(geneDir, "results/")
promoterResDir <- paste0(promoterDir, "results/")
TSSResDir <- paste0(TSSDir, "results/")
TTSResDir <- paste0(TTSDir, "results/")
terminatorResDir <- paste0(terminatorDir, "results/")
genePlotDir <- paste0(geneDir, "histograms/")
promoterPlotDir <- paste0(promoterDir, "histograms/")
TSSPlotDir <- paste0(TSSDir, "histograms/")
TTSPlotDir <- paste0(TTSDir, "histograms/")
terminatorPlotDir <- paste0(terminatorDir, "histograms/")
system(paste0("[ -d ", regionDir, " ] || mkdir ", regionDir))
system(paste0("[ -d ", geneDir, " ] || mkdir ", geneDir))
system(paste0("[ -d ", promoterDir, " ] || mkdir ", promoterDir))
system(paste0("[ -d ", TSSDir, " ] || mkdir ", TSSDir))
system(paste0("[ -d ", TTSDir, " ] || mkdir ", TTSDir))
system(paste0("[ -d ", terminatorDir, " ] || mkdir ", terminatorDir))
system(paste0("[ -d ", geneResDir, " ] || mkdir ", geneResDir))
system(paste0("[ -d ", promoterResDir, " ] || mkdir ", promoterResDir))
system(paste0("[ -d ", TSSResDir, " ] || mkdir ", TSSResDir))
system(paste0("[ -d ", TTSResDir, " ] || mkdir ", TTSResDir))
system(paste0("[ -d ", terminatorResDir, " ] || mkdir ", terminatorResDir))
system(paste0("[ -d ", genePlotDir, " ] || mkdir ", genePlotDir))
system(paste0("[ -d ", promoterPlotDir, " ] || mkdir ", promoterPlotDir))
system(paste0("[ -d ", TSSPlotDir, " ] || mkdir ", TSSPlotDir))
system(paste0("[ -d ", TTSPlotDir, " ] || mkdir ", TTSPlotDir))
system(paste0("[ -d ", terminatorPlotDir, " ] || mkdir ", terminatorPlotDir))

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
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
} else {
  stop("region is not euchromatin, heterochromatin or genomewide")
}

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
  }, mc.cores = detectCores())
  
  # Calculate mean cM/Mb for peak-containing loci (targets) and for
  # each set of peak-less random loci (selected from nontargets)
  targets_cMMbMean <- mean(targets$cMMb, na.rm = T)
  ranLoc_cMMbMean <- unlist(mclapply(seq_along(ranLocGRL), function(x) {
    mean(ranLocGRL[[x]]$cMMb, na.rm = T)
  }, mc.cores = detectCores()))
  
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
