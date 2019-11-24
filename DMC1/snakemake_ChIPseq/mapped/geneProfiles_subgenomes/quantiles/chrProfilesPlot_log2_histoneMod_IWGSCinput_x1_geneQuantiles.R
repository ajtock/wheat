#!/applications/R/R-3.5.0/bin/Rscript

# Plot smoothed library-size-normalized coverage in windows along chromosomes

# Usage:
# ./chrProfilesPlot_log2_histoneMod_IWGSCinput_x1_geneQuantiles.R ASY1_CS ASY1_CS_Rep1_ChIP input H3_input_SRR6350669 both 1Mb 1000000 15 purple4 221119 ASY1_CS_Rep1_ChIP ASY1_CS 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' promoters 4

#markChIPA <- "ASY1_CS"
#libNameChIPA <- "ASY1_CS_Rep1_ChIP"
#markControlA <- "input"
#libNameControlA <- "H3_input_SRR6350669"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#N <- 15
#colourA <- "purple4"
#date <- 221119
#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#region <- "promoters"
#quantiles <- 4

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
markControlA <- args[3]
libNameControlA <- args[4]
align <- args[5]
winName <- args[6]
winSize <- as.numeric(args[7])
N <- as.numeric(args[8])
colourA <- args[9]
date <- args[10]
libName <- args[11]
dirName <- args[12]
featureName <- unlist(strsplit(args[13],
                               split = ","))
region <- args[14]
quantiles <- as.numeric(args[15])

makeTransparent <- function(thisColour, alpha = 210)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}

colourA <- makeTransparent(colourA)

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)
library(plyr)
library(data.table)
library(varhandle)

outDir <- paste0("quantiles_by_log2_", libName,
                 "_control_in_", region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Define plot titles
featureNamePlot <- paste0(toupper(substr(featureName[1], start = 1, stop = 1)),
                          substr(featureName[1], start = 2, stop = 5),
                          " grouped by decreasing ",
                          sub("_\\w+", "", dirName),
                          " in ", region)
# Define quantile names and colours
quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
quantileColours <- c("red", "purple", "blue", "navy")

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
markers <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_mapping_data.txt",
                      header = TRUE)
genes <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA.gff3",
                    colClasses = c(NA,
                                   rep("NULL", 2),
                                   rep(NA, 2),
                                   "NULL", NA, "NULL", NA))
colnames(genes) <- c("chr", "start", "end", "strand", "geneID")
genes <- genes[genes$chr != "chrUn",]
NLRs <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_manually_curated_gene_families/IWGSC_v1.1_nlr_representative_mRNA.gff3",
                   colClasses = c(NA,
                                  rep("NULL", 2),
                                  rep(NA, 2),
                                  "NULL", NA, "NULL", NA))
colnames(NLRs) <- c("chr", "start", "end", "strand", "geneID")
NLRs <- NLRs[NLRs$chr != "chrUn",]

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

## ChIPA
# Calculate log2((ChIP+1)/(Control+1)) coverage within each window
profileChIPAlog2 <- data.frame(chr = as.character(profileChIPA$chr),
                               window = as.numeric(profileChIPA$window),
                               log2CPM = as.numeric(log2((profileChIPA$CPM+1)/(profileControlA$CPM+1))),
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
  filt_chrProfile <- stats::filter(x = chrProfilesChIPA[[x]]$log2CPM,
                                   filter = f,
                                   sides = 2)
  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfilesChIPA[[x]]$chr),
             window = as.integer(chrProfilesChIPA[[x]]$window),
             filt_log2CPM = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfilesChIPA))

minCPM <- min(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    min(filt_chrProfilesChIPA[[x]]$filt_log2CPM)
}, mc.cores = length(filt_chrProfilesChIPA))))
maxCPM <- max(unlist(mclapply(seq_along(filt_chrProfilesChIPA),
  function(x) {
    max(filt_chrProfilesChIPA[[x]]$filt_log2CPM)
}, mc.cores = length(filt_chrProfilesChIPA))))


# Quantiled feature frequency chromosome profiles
featureProfilesList <- mclapply(1:quantiles, function(k) {
  read.table(paste0(outDir,
                    "feature_frequency_per_", winName, "_chromosomeProfiles_",
                    "quantile", k, "_of_", quantiles,
                    "_by_log2_", libName, "_control_in_",
                    region, "_of_",
                    substring(featureName[1][1], first = 1, last = 5), "_in_",
                    paste0(substring(featureName, first = 10, last = 16),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 18), ".tsv"),
             header = T)
})

chrProfilesFeature_quantileList <- mclapply(seq_along(chrs), function(x) {
  lapply(1:quantiles, function(k) {
    featureProfilesList[[k]][featureProfilesList[[k]]$chr == chrs[x],]
  })
}, mc.cores = length(chrs))

filt_chrProfilesFeature_quantileList <- mclapply(seq_along(chrProfilesFeature_quantileList), function(x) {
  lapply(1:quantiles, function(k) {
    filt_chrProfileFeature <- stats::filter(x = chrProfilesFeature_quantileList[[x]][[k]]$winfeatures,
                                            filter = f,
                                            sides = 2)
    filt_chrProfileFeature[1:flank] <- filt_chrProfileFeature[flank+1]
    filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank+1):length(filt_chrProfileFeature)] <- filt_chrProfileFeature[(length(filt_chrProfileFeature)-flank)]
    data.frame(chr = as.character(chrProfilesFeature_quantileList[[x]][[k]]$chr),
               window = as.integer(chrProfilesFeature_quantileList[[x]][[k]]$window),
               filt_feature = as.numeric(filt_chrProfileFeature),
               stringsAsFactors = F)
  })
}, mc.cores = length(chrProfilesFeature_quantileList))

minFeature <- min(unlist(mclapply(seq_along(filt_chrProfilesFeature_quantileList), function(x) {
  min(unlist(lapply(1:quantiles, function(k) {
      min(filt_chrProfilesFeature_quantileList[[x]][[k]]$filt_feature, na.rm = T)
  })))
}, mc.cores = length(filt_chrProfilesFeature_quantileList))))
maxFeature <- max(unlist(mclapply(seq_along(filt_chrProfilesFeature_quantileList), function(x) {
  max(unlist(lapply(1:quantiles, function(k) {
      max(filt_chrProfilesFeature_quantileList[[x]][[k]]$filt_feature, na.rm = T)
  })))
}, mc.cores = length(filt_chrProfilesFeature_quantileList))))


# Plot
pdf(paste0(plotDir, "Wheat_log2_", libNameChIPA, "_", libNameControlA, "_", align, "_",
           "feature_frequency_per_", winName, "_chromosomeProfiles_",
            "quantile", k, "_of_", quantiles,
            "_by_log2_", libName, "_control_in_",
            region, "_of_",
            substring(featureName[1][1], first = 1, last = 5), "_in_",
            paste0(substring(featureName, first = 10, last = 16),
                   collapse = "_"), "_",
            substring(featureName[1][1], first = 18),
           "_chrPlot_winSize", winName, "_smooth", N,
           "_v", date, ".pdf"),
    height = 21, width = 30)
par(mfrow = c(7, 3))
par(mar = c(2.1, 4.5, 2.2, 4.5))
par(mgp = c(3, 1, 0))
for(x in 1:length(filt_chrProfilesChIPA)) {
  chrPlotCov1_featureQuantiles(xplot1 = filt_chrProfilesChIPA[[x]]$window,
                               xplot2 = filt_chrProfilesFeature_quantileList[[x]][[1]]$window,
                               title = paste0(chrs[x], "	", featureNamePlot),
                               cenStart = centromereStart[x],
                               cenEnd = centromereEnd[x],
                               R1End = chrPartitions$R1_R2a[x],
                               R3Start = chrPartitions$R2b_R3[x],
                               rug1 = markers[markers$chromosome == chrs[x],]$physicalPosition,
                               #rug1 = genes[genes$chr == chrs[x],]$start,
                               rug2 = NLRs[NLRs$chr == chrs[x],]$start,
                               regionCol = "red",
                               rug1Col = "grey40",
                               rug2Col = "red3",
                               dat1A = filt_chrProfilesChIPA[[x]]$filt_log2CPM,
                               col1A = colourA,
                               Ylab1 = markChIPA,
                               min1 = -max((minCPM*-1), maxCPM),
                               max1 = max((minCPM*-1), maxCPM),
                               dat2 = filt_chrProfilesFeature_quantileList[[x]],
                               col2 = quantileColours,
                               Ylab2 = "Genes",
                               min2 = minFeature-maxFeature,
                               max2 = maxFeature,
                               legendLoc = "bottomright",
                               legendLabs = quantileNames)
}
dev.off()
