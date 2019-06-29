#!/applications/R/R-3.5.0/bin/Rscript

# Evaluate the relationship between recombination rate and presence of TE classes
# in gene promoters and terminators using permutation tests:
# Divide genes into those whose promoters (1-kb region upstream of TSS)
# or whose terminators (1-kb region downstream of TTS) overlap or do not overlap
# at least one element belonging to a given TE superfamily.
# Those that do not overlap TEs must be within maxDistance of at least one of those that do.
# Calculate mean winName-scaled recombination rate (cM/Mb) for promoters and terminators that
# do and do not (randomSets random sets) overlap TEs and compare to derive significance

# Usage:
#./TE_presence_absence_cMMb_permTest.R 100kb 1 Mariner DTT 500000 10000 0.0001

#winName <- "100kb"
#minMarkerDist <- "1"
#superfamName <- "Mariner"
#superfamCode <- "DTT"
#maxDistance <- 500000
#randomSets <- 10000
#minPval <- 0.0001

args <- commandArgs(trailingOnly = T)
winName <- args[1]
minMarkerDist <- as.numeric(args[2])
superfamName <- args[3]
superfamCode <- args[4]
maxDistance <- as.numeric(args[5])
randomSets <- as.numeric(args[6])
minPval <- as.numeric(args[7])

library(GenomicRanges)
library(parallel)
library(plotrix)

promDir <- "./promoters/"
termDir <- "./terminators/"
promoterResDir <- paste0(promDir, "results/")
terminatorResDir <- paste0(termDir, "results/")
promoterPlotDir <- paste0(promDir, "histograms/")
terminatorPlotDir <- paste0(termDir, "histograms/")
system(paste0("[ -d ", promDir, " ] || mkdir ", promDir))
system(paste0("[ -d ", termDir, " ] || mkdir ", termDir))
system(paste0("[ -d ", promoterResDir, " ] || mkdir ", promoterResDir))
system(paste0("[ -d ", terminatorResDir, " ] || mkdir ", terminatorResDir))
system(paste0("[ -d ", promoterPlotDir, " ] || mkdir ", promoterPlotDir))
system(paste0("[ -d ", terminatorPlotDir, " ] || mkdir ", terminatorPlotDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])

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

# Convert 1000-bp gene promoters and terminators into GRanges
promotersGR <- promoters(genesGR, upstream = 1000, downstream = 0)
promotersGR <- GRanges(promotersGR,
                       geneID = genes$geneID)
print(promotersGR)

source("/projects/ajt200/Rfunctions/TTSplus.R")
terminatorsGR <- TTSplus(genesGR, upstream = -1, downstream = 1000)
terminatorsGR <- GRanges(terminatorsGR,
                         geneID = genes$geneID)
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

# Load table of TEs for given superfamily
# Note that these are in bed format and so the start coordinates are 0-based
TEs <- read.table(paste0(
                  "/home/ajt200/analysis/wheat/featureProfiles/TEs/superfamilies/",
                  "iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
                  superfamName, "_", superfamCode, ".bed"),
                  colClasses = c(rep(NA, 3), rep("NULL", 2), NA))
colnames(TEs) <- c("chr", "start", "end", "strand")
TEs <- TEs[TEs$chr != "chrUn",]
TEsGR <- GRanges(seqnames = TEs$chr,
                 ranges = IRanges(start = TEs$start+1,
                                  end = TEs$end),
                 strand = TEs$strand)

# Obtain gene promoters and terminators that do or don't overlap TEs
# promoters
promoter_TE_overlapsGR <- subsetByOverlaps(x = promotersGR,
                                           ranges = TEsGR,
                                           type = "any",
                                           invert = FALSE,
                                           ignore.strand = TRUE)
promoter_TE_nonoverlapsGR <- subsetByOverlaps(x = promotersGR,
                                              ranges = TEsGR,
                                              type = "any",
                                              invert = TRUE,
                                              ignore.strand = TRUE)
promoter_TE_overlapsGR_extend <- GRanges(seqnames = seqnames(promoter_TE_overlapsGR),
                                         ranges = IRanges(start = start(promoter_TE_overlapsGR)-maxDistance,
                                                          end = end(promoter_TE_overlapsGR)+maxDistance),
                                         strand = strand(promoter_TE_overlapsGR),
                                         cMMb = promoter_TE_overlapsGR$cMMb)
promoter_TE_nonoverlapsGR_nearby <- subsetByOverlaps(x = promoter_TE_nonoverlapsGR,
                                                     ranges = promoter_TE_overlapsGR_extend,
                                                     type = "any",
                                                     invert = FALSE,
                                                     ignore.strand = TRUE)
# terminators
terminator_TE_overlapsGR <- subsetByOverlaps(x = terminatorsGR,
                                             ranges = TEsGR,
                                             type = "any",
                                             invert = FALSE,
                                             ignore.strand = TRUE)
terminator_TE_nonoverlapsGR <- subsetByOverlaps(x = terminatorsGR,
                                                ranges = TEsGR,
                                                type = "any",
                                                invert = TRUE,
                                                ignore.strand = TRUE)
terminator_TE_overlapsGR_extend <- GRanges(seqnames = seqnames(terminator_TE_overlapsGR),
                                           ranges = IRanges(start = start(terminator_TE_overlapsGR)-maxDistance,
                                                            end = end(terminator_TE_overlapsGR)+maxDistance),
                                           strand = strand(terminator_TE_overlapsGR),
                                           cMMb = terminator_TE_overlapsGR$cMMb)
terminator_TE_nonoverlapsGR_nearby <- subsetByOverlaps(x = terminator_TE_nonoverlapsGR,
                                                       ranges = terminator_TE_overlapsGR_extend,
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

# Permutation test function to evaluate if mean cM/Mb at TE-containing loci is
# significantly higher or lower than at random sets of TE-less loci
loci_cMMb_permTest <- function(targets,
                               nontargets,
                               targetsName,
                               targetsNamePlot,
                               resultsDir,
                               plotDir) {
  #targets = promoter_TE_overlapsGR
  #nontargets = promoter_TE_nonoverlapsGR_nearby
  #targetsName = paste0("gene_promoters_overlapping_", superfamName, "_", superfamCode, "_TEs")
  #targetsNamePlot = paste0("Gene promoters overlapping ", superfamName, " (", superfamCode, ") TEs")
  #resultsDir = promoterResDir
  #plotDir = promoterPlotDir
  
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
  
  # Calculate mean cM/Mb for TE-containing loci (targets) and for
  # each set of TE-less random loci (selected from nontargets)
  targets_cMMbMean <- mean(targets$cMMb, na.rm = T)
  ranLoc_cMMbMean <- unlist(mclapply(seq_along(ranLocGRL), function(x) {
    mean(ranLocGRL[[x]]$cMMb, na.rm = T)
  }, mc.cores = detectCores()))
  
  # Determine whether mean cM/Mb values at TE-less random loci are lower than or
  # higher than at TE-containing loci
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
  mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = 1)
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

# Apply to promoters and terminators
loci_cMMb_permTest(targets = promoter_TE_overlapsGR,
                   nontargets = promoter_TE_nonoverlapsGR_nearby,
                   targetsName = paste0("gene_promoters_overlapping_", superfamName, "_", superfamCode, "_TEs"),
                   targetsNamePlot = paste0("Gene promoters overlapping ", superfamName, " (", superfamCode, ") TEs"),
                   resultsDir = promoterResDir,
                   plotDir = promoterPlotDir)
loci_cMMb_permTest(targets = terminator_TE_overlapsGR,
                   nontargets = terminator_TE_nonoverlapsGR_nearby,
                   targetsName = paste0("gene_terminators_overlapping_", superfamName, "_", superfamCode, "_TEs"),
                   targetsNamePlot = paste0("Gene terminators overlapping ", superfamName, " (", superfamCode, ") TEs"),
                   resultsDir = terminatorResDir,
                   plotDir = terminatorPlotDir)
