#!/applications/R/R-3.5.0/bin/Rscript

# Perform hypergeometric tests to determine whether each ASY1 or DMC1
# gene quantile is over-represented or under-represented for
# NLR-encoding genes (as defined by Ksenia Krasileva) and, separately,
# meiotic-module genes (as defined in Alabdullah et al. 2019 Front. Plant Sci.)
# (i.e., is the propotion of NLR-encoding genes or meiotic-module genes
# contained within each gene quantile significantly greater or smaller than expected by chance
# based on the hypergeometric distribution?)

# P-value is the probability of drawing >= length(quantile_NLRs) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_NLRs) [m] + ( length(genome_genes) - length(genome_NLRs)) [n]

### Perform hypergeometric test to determine whether a significant proportion of 
### elements belonging to each transposon family are differentially expressed in kyp suvh5 suvh6

# P-value is the probability of drawing >= length(DE_TEfamIDs) [x] TEs
# in a sample size of length(DE_TEIDs) [k] from a total TE set consisting of
# length(TAIR10_TEfamIDs) [m] + (length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs)) [n]

# Usage 
# ./proportion_NLRs_in_gene_quantiles_hypergeometricTest.R 'ASY1_in_bodies' 'ASY1_CS_Rep1_ChIP' 1 4 'genomewide' 'Agenome_Bgenome_Dgenome' 100000

library(methods)

#quantileBy <- "ASY1_in_bodies"
#libName <- "ASY1_CS_Rep1_ChIP"
#quantileFirst <- 1
#quantileLast <- 4
#region <- "genomewide"
#genomeName <- "Agenome_Bgenome_Dgenome"
#samplesNum <- 100000

args <- commandArgs(trailingOnly = TRUE)
quantileBy <- args[1]
libName <- args[2]
quantileFirst <- as.integer(args[3])
quantileLast <- as.integer(args[4])
region <- args[5]
genomeName <- args[6]
samplesNum <- as.numeric(args[7])

outDir <- paste0("quantiles_by_", quantileBy, "/hypergeometricTests/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load NLRs
NLRs <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/NLRs_Krasileva/NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_representative_mRNA.gff3",
                   header = F, stringsAsFactors = F)
# Replace gene model ID decimal suffix (e.g., ".1")
NLRs$V9 <- sub(pattern = "\\.\\d+", replacement = "",
               x = NLRs$V9)
genome_NLRs <- as.character(NLRs$V9)

#NLRs <- read.table("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/Table_S4_meiotic_GO_genes.tsv",
#                   header = T, stringsAsFactors = F)
#genome_NLRs <- as.character(NLRs$Gene.ID)
#NLRs <- read.table("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/Table_S4_meiotic_gene_orthologs.tsv",
#                   header = T, stringsAsFactors = F)
#genome_NLRs <- as.character(NLRs$Gene.ID)

featuresDF <- read.table(paste0("quantiles_by_", quantileBy,
                                "/WesternEurope/features_", quantileLast, "quantiles_by_",
                                quantileBy, "_of_genes_in_",
                                genomeName, "_", region, "_WesternEurope.txt"),
                         header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
featuresDF$featureID <- sub(pattern = "\\.\\d+", replacement = "",
                            x = featuresDF$featureID)
genome_genes <- featuresDF$featureID
quantile_genes_list <- lapply(quantileFirst:quantileLast, function(x) {
  featuresDF[featuresDF$quantile == paste0("Quantile ", x),]$featureID
})
rm(featuresDF); gc()

# Get the intersection of genome_NLRs and genome_genes
# (this excludes NLRs not assigned to a chromosome)
genome_NLRs <- intersect(genome_NLRs, genome_genes)

# Set class for hypergeometric test results object
setClass("hypergeomTest",
         representation(alternative = "character",
                        alpha0.05 = "numeric",
                        pval = "numeric",
                        observed = "numeric",
                        expected = "numeric",
                        log2obsexp = "numeric",
                        log2alpha = "numeric",
                        quantile_genes = "numeric",
                        proportion_of_quantile = "numeric",
                        random_proportions_of_quantile = "numeric",
                        hypergeomDist = "numeric"))

# P-value is the probability of drawing >= length(quantile_NLRs) [x] features
# in a sample size of length(quantile_genes) [k] from a total feature set consisting of
# length(genome_NLRs) [m] + ( length(genome_genes) - length(genome_NLRs)) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(quantile_NLRs) to length(quantile_genes)

#lapply(seq_along(quantile_genes_list), function(z) {
for(z in seq_along(quantile_genes_list)) {
  quantile_genes <- quantile_genes_list[[z]]
  # Get intersection of gene IDs in quantile z and gene IDs of NLRs
  quantile_NLRs <- intersect(quantile_genes, genome_NLRs)

  # Calculate the P-values for over-representation and under-representation
  # of NLRs among quantile z genes
  set.seed(2847502)
  # Over-representation:
  Pval_overrep <- sum(dhyper(x = length(quantile_NLRs):length(quantile_genes),
                             m = length(genome_NLRs),
                             n = length(genome_genes) - length(genome_NLRs),
                             k = length(quantile_genes)))
  print(Pval_overrep)

  # Or by 1 minus the sum of the probabilities of drawing 0:(length(quantile_NLRs)-1)
  print(1 - sum(dhyper(x = 0:(length(quantile_NLRs)-1),
                       m = length(genome_NLRs),
                       n = length(genome_genes) - length(genome_NLRs),
                       k = length(quantile_genes))))

  # Under-representation
  Pval_underrep <- phyper(q = length(quantile_NLRs),
                          m = length(genome_NLRs),
                          n = length(genome_genes) - length(genome_NLRs),
                          k = length(quantile_genes))
  print(Pval_underrep)

  # Sample without replacement
  hgDist <- rhyper(nn = samplesNum,
                   m = length(genome_NLRs),
                   n = length(genome_genes) - length(genome_NLRs),
                   k = length(quantile_genes))

  # Calculate P-values and significance levels
  if(length(quantile_NLRs) > mean(hgDist)) {
    Pval <- Pval_overrep
    MoreOrLessThanRandom <- "MoreThanRandom"
    alpha0.05 <- quantile(hgDist, probs = 0.95)[[1]]
  } else {
    Pval <- Pval_underrep
    MoreOrLessThanRandom <- "LessThanRandom"
    alpha0.05 <- quantile(hgDist, probs = 0.05)[[1]]
  }

  hgTestResults <- new("hypergeomTest",
                       alternative = MoreOrLessThanRandom,
                       alpha0.05 = alpha0.05,
                       pval = Pval,
                       observed = length(quantile_NLRs),
                       expected = mean(hgDist),
                       log2obsexp = log2( length(quantile_NLRs) / mean(hgDist) ),
                       log2alpha  = log2( alpha0.05 / mean(hgDist) ),
                       quantile_genes = length(quantile_genes),
                       proportion_of_quantile = length(quantile_NLRs) / length(quantile_genes),
                       random_proportions_of_quantile = hgDist / length(quantile_genes),
                       hypergeomDist = hgDist)
  save(hgTestResults,
       file = paste0(outDir,
                     orderingFactor, "_", pop_name[x], "_random_nonIDs_permTest_for_", gsub(" ", "_", featureNamePl
ot),
                     "_in_quantile", quantileNo, "_of_", quantiles,
                     "_by_log2_", libName, "_control_in_", region, "_of_",
                     substring(featureName[1][1], first = 1, last = 5), "_in_",
                     paste0(substring(featureName, first = 10, last = 16),
                            collapse = "_"), "_",
                     substring(featureName[1][1], first = 18),
                     ".RData"))


  save(hgTestResults,
       file = paste0(outDir,
                     baseName,
                     "_hypergeometricTestResults_",
                     TEfamNames[x],
                     ".RData"))

  library(plotrix)
  # Generate histogram
  pdf(paste0(plotDir,
             "hist_",
             baseName,
             "_hypergeometricTestResults_",
             TEfamNames[x],
             ".pdf"),
             height = 4, width = 5)
  ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
  #options(scipen = 100)
  # Calculate max density
  maxDensityPlus <- max(density(hgTestResults@random_proportions_DE_TEfam)$y)*1.2
  if(hgTestResults@proportion_DE_TEfam > mean(hgTestResults@random_proportions_DE_TEfam)) {
    alpha0.05 <- quantile(hgTestResults@random_proportions_DE_TEfam, 0.95)[[1]]
    xlim <- c(pmax(0, min(hgTestResults@random_proportions_DE_TEfam)-.1),
              pmax(hgTestResults@proportion_DE_TEfam+.1, alpha0.05+.1))
    Pval <- Pval_overrep
  } else {
    alpha0.05 <- quantile(hgTestResults@random_proportions_DE_TEfam, 0.05)[[1]]
    xlim <- c(pmax(0, pmin(hgTestResults@proportion_DE_TEfam-.1, alpha0.05-.1)),
              max(hgTestResults@random_proportions_DE_TEfam)+.1)
    Pval <- Pval_underrep
  }
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  hist(hgTestResults@random_proportions_DE_TEfam,
       freq = FALSE,
       col = "grey70",
       border = NA,
       lwd = 2,
       xlim = xlim,
       ylim = c(0,
                maxDensityPlus),
       xlab = paste0("Proportion of ", regulated, " transposable elements \n that are ",
                     TEfamNamesPlot[x], " elements"),
       ylab = "Density",
       main = "",
       cex.lab = 1, cex.axis = 1)
  titleText <- list(bquote(.(baseName)),
                    bquote(italic("P")~" = "~.(Pval)),
                    bquote("Samples (hypergeometric distribution) = "~.(samplesChar)))
  mtext(do.call(expression, titleText), side = 3, line = 2:0, cex = 1)
  lines(density(hgTestResults@random_proportions_DE_TEfam), lwd = 1.5)
  ablineclip(v = mean(hgTestResults@random_proportions_DE_TEfam),
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
  ablineclip(v = hgTestResults@proportion_DE_TEfam,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
  ablineclip(v = alpha0.05,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
  text(x = c(pmax(0.05, min(hgTestResults@random_proportions_DE_TEfam)-.05),
             mean(hgTestResults@random_proportions_DE_TEfam),
             hgTestResults@proportion_DE_TEfam,
             alpha0.05),
       y = c(maxDensityPlus*.95,
             maxDensityPlus,
             maxDensityPlus,
             maxDensityPlus*.95),
       labels = c("Simulated",
                  "Expected",
                  "Observed",
                  expression(alpha~" = 0.05")),
       col = c("grey70",
               "black",
               "forestgreen",
               "red"),
       cex = 0.7)
  box(lwd = 2)
  dev.off()

})
