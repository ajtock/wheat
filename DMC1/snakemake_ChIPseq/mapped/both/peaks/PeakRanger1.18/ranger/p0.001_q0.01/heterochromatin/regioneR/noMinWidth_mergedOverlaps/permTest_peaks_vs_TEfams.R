#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 100G -c 48 "/applications/R/R-3.3.2/bin/Rscript permTest_peaks_vs_TEfams.R DMC1_Rep1_ChIP 10000 'heterochromatin' 'A'"

args <- commandArgs(trailingOnly = TRUE)
libNameChIP <- args[1]
perms <- as.numeric(args[2])
region <- args[3]
genomeName <- args[4]

library(regioneR)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = chrStart,
                                     end = chrLens),
                    strand = "*")
genomeGR <- genomeGR[grep(genomeName,
                          seqnames(genomeGR))@values]

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
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin or genomewide")
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
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin or genomewide")
}

# Import peaks as GRanges object
load(paste0("../../../", libNameChIP,
            "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
peaksGR <- peaksGR[grep(genomeName,
                        seqnames(peaksGR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaksGR <- peaksGR[-subjectHits(mask_peaks_overlap)]
strand(peaksGR) <- "*"
print("***********peaks***********")
print(peaksGR)

# TE superfamilies
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
  TEsGRtmp <- TEsGRtmp[grep(genomeName,
                            seqnames(TEsGRtmp))@values]
  # Subset to include only those not overlapping masked region
  mask_TEstmp_overlap <- findOverlaps(query = maskGR,
                                      subject = TEsGRtmp,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
  TEsGRtmp[-subjectHits(mask_TEstmp_overlap)]
})

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(845934)
ptPeaksOtherPerChrom <- lapply(seq_along(othersGRL), function(x) {
  permTest(A = peaksGR,
           B = othersGRL[[x]],
           genome = genomeGR,
           mask = maskGR,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE,
           per.chromosome = TRUE,
           evaluate.function = numOverlaps,
           count.once = TRUE,
           ntimes = perms,
           mc.set.seed = FALSE,
           mc.cores = detectCores())
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0("permTest_", as.character(perms), "perms_",
                   libNameChIP, "_peaks_vs_TEfams_in_",
                   genomeName, "genome_", region,
                   ".RData"))

# Summarise results in a table
featureName <- NULL
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksOtherPerChrom)) {
  featureNamei <- print(otherNames[i])
  featureName <- c(featureName, featureNamei)
  noOfFeaturesi <- print(length(othersGRL[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksOtherPerChromDataFrame <- cbind(featureName, noOfFeatures, expected, observed, pval, zscore)
colnames(ptPeaksOtherPerChromDataFrame) <- c("feature", "n", "expected", "observed", "pval", "zscore")
write.table(ptPeaksOtherPerChromDataFrame,
            file = paste0("permTest_", as.character(perms), "perms_",
                          libNameChIP, "_peaks_vs_TEfams_in_",
                          genomeName, "genome_", region,
                          "_DataFrame.txt"),
            sep = "\t", quote = F, row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(paste0(plotDir, otherNames[i],
             "_permTest_", as.character(perms), "perms_",
             libNameChIP, "_peaks_in_",
             genomeName, "genome_", region,
             "_perChrom.pdf"),
      width = 10, height = 7)
  plot(ptPeaksOtherPerChrom[[i]], main = paste0(libNameChIP, " peaks vs ", otherNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = othersGRL[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = othersGRL[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = othersGRL[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(paste0(plotDir, otherNames[i],
             "_localZscore_permTest_", as.character(perms), "perms_",
             libNameChIP, "_peaks_in_",
             genomeName, "genome_", region,
             "_w1kb_s50bp_w10kb_s500bp_w",
             win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(libNameChIP, " peaks vs ", otherNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(libNameChIP, " peaks vs ", otherNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(libNameChIP, " peaks vs ", otherNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(libNameChIP, " peaks vs ", otherNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(libNameChIP, " peaks vs ", otherNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(libNameChIP, " peaks vs ", otherNames[i], " (~", win, "-kb shift)"))
  dev.off()
}
