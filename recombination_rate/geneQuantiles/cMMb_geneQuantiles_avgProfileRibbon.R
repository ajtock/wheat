#!/applications/R/R-3.5.0/bin/Rscript

# Divide genes into quantiles based on mean winName-scaled
# recombination rate (cM/Mb) from promoters to terminators.
# Plot average coverage profiles for a given dataset
# around these genes separated into cM/Mb quantiles.

# Usage:
#./cMMb_geneQuantiles_avgProfileRibbon.R 100kb 1 4 ASY1_CS ASY_CS_Rep1_ChIP input H3_input_SRR6350669 both purple4 genes 20bp 20 2kb 2000 'euchromatin'

winName <- "100kb"
minMarkerDist <- 1
quantiles <- 4
markChIP <- "ASY1_CS"
libNameChIP <- "ASY1_CS_Rep1_ChIP"
markControl <- "input"
libNameControl <- "H3_input_SRR6350669"
align <- "both"
colour <- "purple4"
featureName <- "genes"
binName <- "20bp"
binSize <- 20
flankName <- "2kb"
flankSize <- 2000
region <- "euchromatin"

args <- commandArgs(trailingOnly = T)
winName <- args[1]
minMarkerDist <- as.numeric(args[2])
quantiles <- as.numeric(args[3])
markChIP <- args[4]
libNameChIP <- args[5]
markControl <- args[6]
libNameControl <- args[7]
align <- args[8]
colour <- args[9]
featureName <- args[10]
binName <- args[11]
binSize <- as.numeric(args[12])
flankName <- args[13]
flankSize <- as.numeric(args[14])
region <- args[15]

library(GenomicRanges)
library(parallel)

regionDir <- paste0(region, "/")
plotDir <- paste0(regionDir, "plots/")
system(paste0("[ -d ", regionDir, " ] || mkdir ", regionDir))
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

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
} else if(region == "genomewide") {
  maskGR <- GRanges()
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
# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = maskGR,
                                   subject = genesGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
genesGR <- genesGR[-subjectHits(mask_genes_overlap)]
# Extend gene boundaries to include promoters and terminators for calculation of
# winName-scaled recombination rate
genesGR <- GRanges(seqnames = seqnames(genesGR),
                   ranges = IRanges(start = start(genesGR)-1000,
                                    end = end(genesGR)+1000),
                   strand = strand(genesGR),
                   geneID = genesGR$geneID)
print(genesGR)

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

# Obtain winName-scaled cMMb values for each gene between promoter and terminator
# Where features overlap more than one winName window, calculate mean cMMb
gene_cMMb_overlaps <- findOverlaps(query = genesGR,
                                   subject = cMMbGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
gene_cMMb_overlapsList <- lapply(seq_along(genesGR), function(x) {
  subjectHits(gene_cMMb_overlaps)[queryHits(gene_cMMb_overlaps) == x]
})
## OR
#gene_cMMb_overlapsList <- getOverlaps(coordinates = genesGR,
#                                      segments = cMMbGR,
#                                      overlapType = "overlapping",
#                                      whichOverlaps = TRUE,
#                                      ignoreStrand = TRUE)
gene_cMMb <- sapply(gene_cMMb_overlapsList,
                    function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
genesGR <- GRanges(genesGR,
                   geneID = genesGR$geneID,
                   cMMb = gene_cMMb)

# Separate genes into genomes
genomeNames <- c("A", "B", "D", "")
genomeNamesPlot <- c("A", "B", "D", "all")
genesGR_genomeList <- lapply(seq_along(genomeNames), function(z) {
  genesGR[grep(genomeNames[z],
               seqnames(genesGR))]
})

# Order genes by decreasing cM/Mb
# and divide into quantiles for average coverage profiling
sapply(seq_along(genesGR_genomeList), function(z) {
  genesDF <- data.frame(genesGR_genomeList[[z]])
  genesDF_ordered <- genesDF[order(genesDF$cMMb,
                                   decreasing = T,
                                   na.last = NA),]
  quantilesN <- round(dim(genesDF_ordered)[1]/quantiles)
  quantilesCum <- cumsum(c(1,
                           rep(quantilesN,
                               times = quantiles)))
  if(quantilesCum[length(quantilesCum)] < dim(genesDF_ordered)[1]) {
    quantilesCum <- c(quantilesCum, dim(genesDF_ordered)[1])
  }
  quantilesStats <- data.frame()
  for(j in 1:length(quantilesCum)) {
    print(j)
    if(j == 1) {
      quantilejFeatures <- genesDF_ordered[(quantilesCum[j]):(quantilesCum[j+1]-1),]
      print(paste0("condition 1: quantile ", j))
      print(dim(quantilejFeatures))
    }
    if(j > 1 & j < quantiles) {
      quantilejFeatures <- genesDF_ordered[(quantilesCum[j]):(quantilesCum[j+1]-1),]
      print(paste0("condition 2: quantile ", j))
      print(dim(quantilejFeatures))
    }
    if(j == quantiles) {
      quantilejFeatures <- genesDF_ordered[(quantilesCum[j]):(quantilesCum[length(quantilesCum)]),]
      print(paste0("condition 3: quantile ", j))
      print(dim(quantilejFeatures))
    }
    if(j <= quantiles) {
      stats <-  data.frame(quantile = as.integer(j),
                           n = as.integer(dim(quantilejFeatures)[1]),
                           mean_width = as.integer(round(mean(quantilejFeatures$width))),
                           total_width = as.integer(sum(quantilejFeatures$width)),
                           mean_cMMb = as.numeric(mean(quantilejFeatures$cMMb)))
    quantilesStats <- rbind(quantilesStats, stats)
    write.table(quantilejFeatures,
                file = paste0(regionDir,
                              "genes_in_", genomeNames[z], "genome_", region,
                              "_quantile", j, "_of_", quantiles, "quantiles_", 
                              "ordered_by_", winName, "Scaled_cMMb_",
                              "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    }
  }
  write.table(quantilesStats,
              file = paste0(regionDir,
                            "summary_genes_in_", genomeNames[z], "genome_", region, "_",
                            quantiles, "quantiles_", 
                            "ordered_by_", winName, "Scaled_cMMb_",
                            "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
})

genesDF_genomeList_quantiles <- lapply(seq_along(genesGR_genomeList),
  function(z) {
    lapply(seq_along(1:quantiles), function(j) {
      read.table(paste0(regionDir,
                        "genes_in_", genomeNames[z], "genome_", region,
                        "_quantile", j, "_of_", quantiles, "quantiles_", 
                        "ordered_by_", winName, "Scaled_cMMb_",
                        "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                 header = TRUE)
    })
})



# Load feature coverage matrices
# ChIP
if(libNameChIP %in% c("H3K4me3_ChIP_SRR6350668",
                      "H3K27me3_ChIP_SRR6350666",
                      "H3K36me3_ChIP_SRR6350670",
                      "H3K9ac_ChIP_SRR6350667",
                      "CENH3_ChIP_SRR1686799")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                       markChIP, "/snakemake_ChIPseq/mapped/geneProfiles/matrices/")
} else {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/",
                       markChIP, "/snakemake_ChIPseq/mapped/geneProfiles/matrices/")
}
covMatChIP <- read.table(paste0(covDirChIP,
                                libNameChIP,
                                "_MappedOn_wheat_v1.0_lowXM_",
                                align, "_sort_norm_",
                                featureName, "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         skip = 3)
covMatChIP <- data.frame(geneID = genes$geneID,
                         covMatChIP)
covMatChIP_genomeList_quantiles <- lapply(seq_along(genesDF_genomeList_quantiles),
  function(z) {
    lapply(seq_along(1:quantiles), function(j) {
      covMatChIP[covMatChIP$geneID %in% genesDF_genomeList_quantiles[[z]][[j]]$geneID,]
    })
})

# Control
if(libNameControl == "MNase_Rep1") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/",
                          "MNase/snakemake_ChIPseq/mapped/", align, "/bg/")
  profileControl <- read.table(paste0(covDirControl, "MNase_Rep1_MappedOn_wheat_v1.0_lowXM_",
                                      align, "_sort_norm_binSize", winName, ".bedgraph"))
} else if(libNameControl == "H3_input_SRR6350669") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                          "input/snakemake_ChIPseq/mapped/", align, "/bg/")
  profileControl <- read.table(paste0(covDirControl, "H3_input_SRR6350669_MappedOn_wheat_v1.0_lowXM_",
                                      align, "_sort_norm_binSize", winName, ".bedgraph"))
} else {
  if(!(libNameControl %in% c("MNase_Rep1", "H3_input_SRR6350669"))) {
    stop("libNameControl is neither MNase_Rep1 nor H3_input_SRR6350669")
  }
}

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
                  expression(alpha~" = 0.05")),
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
