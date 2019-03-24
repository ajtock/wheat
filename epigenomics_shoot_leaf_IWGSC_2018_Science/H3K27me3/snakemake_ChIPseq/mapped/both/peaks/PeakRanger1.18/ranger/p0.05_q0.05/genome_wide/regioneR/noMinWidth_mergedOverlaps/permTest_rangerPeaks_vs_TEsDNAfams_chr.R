#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 100G -c 48 "/applications/R/R-3.3.2/bin/Rscript permTest_rangerPeaks_vs_TEsDNAfams_chr.R H3K4me3_Rep1_ChIP chr1A"

args <- commandArgs(trailingOnly = TRUE)
ChIPLibName <- args[1]
chrName <- args[2]

library(regioneR)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,3])
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
genome <- genome[seqnames(genome) == chrName]

system(paste0("[ -d ", chrName, " ] || mkdir ", chrName))
system(paste0("[ -d ", chrName, "/plots ] || mkdir ", chrName, "/plots"))
system(paste0("[ -d ", chrName, "/plots/TEsDNA ] || mkdir ", chrName, "/plots/TEsDNA")
)
outDir <- paste0("./", chrName, "/")
DNAplotDir <- paste0(outDir, "plots/TEsDNA/")

# Import peaks as GRanges object
load(paste0("../../../", ChIPLibName,
            "_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData"))
peaksGR <- rangerPeaksGRmergedOverlaps
strand(peaksGR) <- "*"
peaksGR <- peaksGR[seqnames(peaksGR) == chrName]
print("***********peaks***********")
print(peaksGR)
print(length(peaksGR))

DNAfamNames <-  c(
                  "CACTA_DTC",
                  "Harbinger_DTH",
                  "hAT_DTA",
                  "Helitrons_DHH",
                  "Mariner_DTT",
                  "Mutator_DTM",
                  "Unclassified_class_2_DXX",
                  "Unclassified_with_TIRs_DTX"
                 )
RNAfamNames <-  c(
                  "Copia_LTR_RLC",
                  "Gypsy_LTR_RLG",
                  "LINE_RIX",
                  "SINE_SIX",
                  "Unclassified_LTR_RLX",
                  "Unclassified_repeats_XXX"
                 )
TEdir <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/superfamilies/"

### DNA TEs
TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(paste0(TEdir,
                              "iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
                              DNAfamNames[x], ".bed"), header = F)
  TEsDNAGRtmp <- GRanges(seqnames = TEsDNA$V1,
                         ranges = IRanges(start = TEsDNA$V2+1,
                                          end = TEsDNA$V3),
                         strand = "*")
  TEsDNAGRtmp[seqnames(TEsDNAGRtmp) == chrName]
})

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(845934)
ptPeaksTEsDNAPerChrom <- lapply(seq_along(TEsDNAGR), function(x) {
  permTest(A = peaksGR, B = TEsDNAGR[[x]], genome = genome,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = detectCores())
})

for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
  assign(paste0(DNAfamNames[i]), ptPeaksTEsDNAPerChrom[[i]])
}
save(ptPeaksTEsDNAPerChrom,
     file = paste0(outDir,
                   "permTest_", ChIPLibName,
                   "_rangerPeaks_vs_TEsDNA_",
                   chrName, ".RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
  noOfFeaturesi <- print(length(TEsDNAGR[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksTEsDNAPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptPeaksTEsDNAPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_", ChIPLibName,
                          "_rangerPeaks_vs_TEsDNA_",
                          chrName, "_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
  pdf(paste0(DNAplotDir, DNAfamNames[i],
             "_permTest_nperm10000_", ChIPLibName,
             "_rangerPeaks_perChrom_", chrName, ".pdf"),
      width = 10, height = 7)
  plot(ptPeaksTEsDNAPerChrom[[i]], main = paste0(ChIPLibName, " rangerPeaks vs ", DNAfamNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGR[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGR[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGR[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(paste0(DNAplotDir, DNAfamNames[i],
             "_localZscore_permTest_nperm10000_"
             ChIPLibName, "_rangerPeaks_w1kb_s50bp_w10kb_s500bp_w",
             win ,"kb_s", step, "bp_perChrom_", chrName, ".pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(ChIPLibName, " rangerPeaks vs ", DNAfamNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(ChIPLibName, " rangerPeaks vs ", DNAfamNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(ChIPLibName, " rangerPeaks vs ", DNAfamNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(ChIPLibName, " rangerPeaks vs ", DNAfamNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(ChIPLibName, " rangerPeaks vs ", DNAfamNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(ChIPLibName, " rangerPeaks vs ", DNAfamNames[i], " (~", win, "-kb shift)"))
  dev.off()
}
