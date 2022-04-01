#!/usr/bin/env Rscript

# Usage:
# ./AxC_mapped_marker_interval_correlation_matrix.R

plotDir <- paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

library(GenomicRanges)
library(Hmisc) # includes rcorr() function which computes significance levels for Pearson and Spearman correlations
library(reshape)
library(ggplot2)
library(ggcorrplot)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

# Load cM/Mb data
dat <- read.table("../AxC_mapped_marker_intervals_mean_ChIPseq_DNAmethyl_cMMb.tsv", header = T)
colnames(dat)
# [1] "chr"                      "start"
# [3] "end"                      "width"
# [5] "physical_marker"          "wt_marker"
# [7] "wt_cM"                    "fancm_marker"
# [9] "fancm_cM"                 "log2_ASY1_CS_input"
#[11] "log2_DMC1_input"          "log2_H3K4me1_input"
#[13] "log2_H3K4me3_MNase"       "log2_H3K27ac_input"
#[15] "log2_H3K27me3_input"      "log2_H3K36me3_input"
#[17] "log2_H3K9me2_MNase"       "log2_H3K27me1_MNase"
#[19] "log2_CENH3_input"         "mCpG"
#[21] "mCHG"                     "mCHH"
#[23] "cMMb_1Mb"                 "cMMb_10Mb"
#[25] "wt_cM_inter"              "fancm_cM_inter"
#[27] "wt_cMMb_inter"            "fancm_cMMb_inter"
#[29] "diff_fancm_wt_cM_inter"   "diff_fancm_wt_cMMb_inter"
#[31] "divi_fancm_wt_cM_inter"   "divi_fancm_wt_cMMb_inter"
#[33] "l2fc_fancm_wt_cM_inter"   "l2fc_fancm_wt_cMMb_inter"

colnames(dat) <- c(colnames(dat)[1:9],
                   "ASY1", "DMC1", "H3K4me1", "H3K4me3", "H3K27ac",
                   "H3K27me3", "H3K36me3", "H3K9me2", "H3K27me1",
                   "CENH3", "mCG", "mCHG", "mCHH", "cMMb_1Mb", "cMMb_10Mb",
                   colnames(dat)[25:ncol(dat)])

profilesDF <- cbind(dat$diff_fancm_wt_cMMb_inter,
                    dat[,10:24])
profilesDF <- profilesDF[,-15]
profilesDF <- data.frame("Diff_cMMb" = profilesDF[,1],
                         "IWGSC_cMMb" = profilesDF[,ncol(profilesDF)],
                         profilesDF[,c(-1,-ncol(profilesDF))])

profileNames <- colnames(profilesDF)

# Create correlation matrix
corMat <- round(cor(profilesDF,
                    method = "spearman",
                    use = "pairwise.complete.obs"),
                digits = 2)

# Set duplicates to NA
for(x in 1:dim(corMat)[1]) {
  corMat[x, x] <- NA
  if(x > 1) {
    corMat[x, 1:x-1] <- NA
  }
}
corMat <- corMat[,-1]

# Convert into reshape::melt formatted data.frame
# and remove duplicate pairs
corDat <- melt(corMat)
corDat <- corDat[-which(is.na(corDat[,3])),]

# Order the data.frame for plotting
profileNamesList <- as.list(profileNames)
names(profileNamesList) <- profileNames
levels(corDat$X1) <- rev(profileNamesList) 
levels(corDat$X2) <- profileNamesList[-1]

# Get P-values for correlation matrix
corMatSig <- rcorr(as.matrix(profilesDF),
                   type = "spearman")$P
# Set duplicates to NA
for(x in 1:dim(corMatSig)[1]) {
  corMatSig[x, x] <- NA
  if(x > 1) {
    corMatSig[x, 1:x-1] <- NA
  }
}
corMatSig <- corMatSig[,-1]

# Convert into reshape::melt formatted data.frame
# and remove duplicate pairs
corDatSig <- melt(corMatSig)
corDatSig <- corDatSig[-which(is.na(corDatSig[,3])),]

# Standardise P-values to a sample size of 100 (q-values) as proposed by
# Good (1982) Standardized tail-area probabilities. Journal of Computation and Simulation 16: 65-66
# and summarised by Woolley (2003):
# https://stats.stackexchange.com/questions/22233/how-to-choose-significance-level-for-a-large-data-set
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.518.5341&rep=rep1&type=pdf
# Woolley (2003): "Clearly, the meaningfulness of the p-value diminishes as the sample size increases";
# Anne Z. (2012, Pearson eCollege, Denver): "In the real world, there are unlikely to be semi-partial correlations
# that are exactly zero, which is the null hypothesis in testing significance of a regression coefficient."
# Formally, the standardised p-value is defined as:
# q = min(0.5, p * sqrt( (n/100) ))
# Woolley (2003): "The value of 0.5 is somewhat arbitrary, though its purpose is to avoid q-values of greater than 1."
n <- dim(profilesDF)[1]
corDatSig$value <- sapply(corDatSig$value, function(x) {
  round(min(0.5, x * sqrt( (n/100) )),
        digits = 2)
})

#corDatSig$value <- sapply(corDatSig$value, function(x) {
#  round(x, digits = 2)
#})

# Order the data.frame for plotting
levels(corDatSig$X1) <- rev(profileNamesList) 
levels(corDatSig$X2) <- profileNamesList[-1]

# Plot
ggObj <- ggplot(data = corDat,
                mapping = aes(X2, X1, fill = value)) +
  geom_tile() +
#  geom_text(mapping = aes(X2, X1, label = value), size = 5) +
  geom_text(data = corDatSig,
            mapping = aes(X2, X1, label = value), size = 8) +
  scale_fill_gradient2(name = bquote("Spearman's" ~ italic(r[s])),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, breaks = seq(-1, 1, by = 0.4), limits = c(-1, 1)) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(barwidth = 40, barheight = 6,
                                title.position = "top", title.hjust = 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0, size = 39, colour = "black"),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 39, colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.justification = c(1, 0),
        legend.position = c(0.65, 0.05),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(5.5, 70.5, 5.5, 5.5), "pt"),
        plot.title = element_text(hjust = 0.5, size = 30, colour = "black")) +
  ggtitle(bquote("Spearman's" ~ italic(r[s]) ~ "for AxC mapped marker intervals"))
ggsave(paste0(plotDir,
              "Spearman_correlation_matrix_AxC_mapped_marker_intervals_mean_ChIPseq_DNAmethyl_cMMb_qVals.pdf"),
       plot = ggObj, height = 20, width = 20)
