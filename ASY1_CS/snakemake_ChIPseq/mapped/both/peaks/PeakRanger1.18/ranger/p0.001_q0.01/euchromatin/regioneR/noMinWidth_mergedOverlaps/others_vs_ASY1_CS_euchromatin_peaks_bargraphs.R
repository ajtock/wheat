#!/applications/R/R-3.5.0/bin/Rscript

# Plot bar graph of log2(observed:expected) peaks overlapping other features

# Usage:
# ./others_vs_ASY1_CS_euchromatin_peaks_bargraphs.R "Euchromatin ASY1 peaks" "A genome" "B genome" "D genome" Agenome Bgenome Dgenome Agenome euchromatin 10000

library(ggplot2)
library(ggthemes)

#plotTitle <- "Euchromatin ASY1 peaks"
#pt1Name <- "A genome"
#pt2Name <- "B genome"
#pt3Name <- "D genome"
#pt1LibName <- "Agenome" 
#pt2LibName <- "Bgenome"
#pt3LibName <- "Dgenome"
#ptOrderLibName <- "Agenome"
#region <- "euchromatin" 
## Number of permutations (randomisations) performed
#perms <- "10000"

args <- commandArgs(trailingOnly = T)
plotTitle <- args[1]
pt1Name <- args[2]
pt2Name <- args[3]
pt3Name <- args[4]
pt1LibName <- args[5]
pt2LibName <- args[6]
pt3LibName <- args[7]
ptOrderLibName <- args[8]
region <- args[9]
# Number of permutations (randomisations) performed
perms <- as.character(args[10])

ptOrderDir <- paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                     "euchromatin/regioneR/noMinWidth_mergedOverlaps/")
pt1Dir <- paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                 region, "/regioneR/noMinWidth_mergedOverlaps/")
pt2Dir <- paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                 region, "/regioneR/noMinWidth_mergedOverlaps/")
pt3Dir <- paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                 region, "/regioneR/noMinWidth_mergedOverlaps/")

plotDir <- "./bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

otherNamesPlot <- c(
                    "H3K4me3",
                    "H3K9me2",
                    "H3K27me1",
                    "H3K27me3",
                    "H3K36me3",
                    "H3K9ac",
                    "H2A.Z",
                    "Nucleosomes",
                    "Genes",
                    "Gene promoters",
                    "Gene 5' ends",
                    "Gene 3' ends",
                    "Gene terminators",
                    "NLRs",
                    "NLR promoters",
                    "NLR 5' ends",
                    "NLR 3' ends",
                    "NLR terminators"
                   )

# Load permutation test results for peak set to be used for ordering
# of other features in bar graph
load(paste0(ptOrderDir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_others_in_",
            ptOrderLibName, "_euchromatin.RData"))
ptOrder <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

# Load permutation test results to be used for plotting
load(paste0(pt1Dir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_others_in_",
            pt1LibName, "_", region, ".RData"))
pt1 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

load(paste0(pt2Dir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_others_in_",
            pt2LibName, "_", region, ".RData"))
pt2 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

load(paste0(pt3Dir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_others_in_",
            pt3LibName, "_", region, ".RData"))
pt3 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

# ptOrder
ptOrder_Pval <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$pval
})
ptOrder_Obs <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$observed
})
ptOrder_Perm <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$permuted
})
ptOrder_Exp <- lapply(seq_along(ptOrder), function(x) {
  mean(ptOrder[[x]]$numOverlaps$permuted)
})
ptOrder_log2ObsExp <- lapply(seq_along(ptOrder_Obs), function(x) {
  log2(ptOrder_Obs[[x]]/ptOrder_Exp[[x]])
})
ptOrder_Zscore <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$zscore
})
ptOrder_AltHyp <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$alternative
})
ptOrder_alpha0.05 <- lapply(seq_along(ptOrder_Perm), function(x) {
  if(ptOrder_AltHyp[[x]] == "greater") {
    quantile(ptOrder_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(ptOrder_Perm[[x]], 0.05)[[1]]
  }
})
ptOrder_log2alpha0.05 <- lapply(seq_along(ptOrder_alpha0.05), function(x) {
  log2(ptOrder_alpha0.05[[x]]/ptOrder_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
ptOrder_log2ObsExp_sorted <- sort.int(unlist(ptOrder_log2ObsExp),
                                      decreasing = T)
ptOrder_log2alpha0.05_sorted <- unlist(ptOrder_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                                      decreasing = T,
                                                                      index.return = T)$ix])
#ptOrder_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                             decreasing = T,
#                                             index.return = T)$ix]
ptOrder_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                     decreasing = T,
                                                     index.return = T)$ix]

# pt1
pt1_Pval <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$pval
})
pt1_Obs <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$observed
})
pt1_Perm <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$permuted
})
pt1_Exp <- lapply(seq_along(pt1), function(x) {
  mean(pt1[[x]]$numOverlaps$permuted)
})
pt1_log2ObsExp <- lapply(seq_along(pt1_Obs), function(x) {
  log2(pt1_Obs[[x]]/pt1_Exp[[x]])
})
pt1_Zscore <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$zscore
})
pt1_AltHyp <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$alternative
})
pt1_alpha0.05 <- lapply(seq_along(pt1_Perm), function(x) {
  if(pt1_AltHyp[[x]] == "greater") {
    quantile(pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_Perm[[x]], 0.05)[[1]]
  }
})
pt1_log2alpha0.05 <- lapply(seq_along(pt1_alpha0.05), function(x) {
  log2(pt1_alpha0.05[[x]]/pt1_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt1_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt1_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# pt2
pt2_Pval <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$pval
})
pt2_Obs <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$observed
})
pt2_Perm <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$permuted
})
pt2_Exp <- lapply(seq_along(pt2), function(x) {
  mean(pt2[[x]]$numOverlaps$permuted)
})
pt2_log2ObsExp <- lapply(seq_along(pt2_Obs), function(x) {
  log2(pt2_Obs[[x]]/pt2_Exp[[x]])
})
pt2_Zscore <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$zscore
})
pt2_AltHyp <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$alternative
})
pt2_alpha0.05 <- lapply(seq_along(pt2_Perm), function(x) {
  if(pt2_AltHyp[[x]] == "greater") {
    quantile(pt2_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2_Perm[[x]], 0.05)[[1]]
  }
})
pt2_log2alpha0.05 <- lapply(seq_along(pt2_alpha0.05), function(x) {
  log2(pt2_alpha0.05[[x]]/pt2_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt2_log2ObsExp_sorted <- unlist(pt2_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt2_log2alpha0.05_sorted <- unlist(pt2_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt2_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt2_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]
# pt3
pt3_Pval <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$pval
})
pt3_Obs <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$observed
})
pt3_Perm <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$permuted
})
pt3_Exp <- lapply(seq_along(pt3), function(x) {
  mean(pt3[[x]]$numOverlaps$permuted)
})
pt3_log2ObsExp <- lapply(seq_along(pt3_Obs), function(x) {
  log2(pt3_Obs[[x]]/pt3_Exp[[x]])
})
pt3_Zscore <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$zscore
})
pt3_AltHyp <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$alternative
})
pt3_alpha0.05 <- lapply(seq_along(pt3_Perm), function(x) {
  if(pt3_AltHyp[[x]] == "greater") {
    quantile(pt3_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt3_Perm[[x]], 0.05)[[1]]
  }
})
pt3_log2alpha0.05 <- lapply(seq_along(pt3_alpha0.05), function(x) {
  log2(pt3_alpha0.05[[x]]/pt3_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt3_log2ObsExp_sorted <- unlist(pt3_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt3_log2alpha0.05_sorted <- unlist(pt3_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt3_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt3_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# Combine in data.frame
df <- data.frame(Sample = rep(c(pt1Name, pt2Name, pt3Name),
                              each = length(ptOrder_log2ObsExp_sorted)),
                 Feature = rep(pt1_otherNamesPlot_sorted, 3),
                 log2ObsExp = c(pt1_log2ObsExp_sorted,
                                pt2_log2ObsExp_sorted,
                                pt3_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_log2alpha0.05_sorted,
                                   pt2_log2alpha0.05_sorted,
                                   pt3_log2alpha0.05_sorted))

df$Feature <- factor(df$Feature,
                               levels = pt1_otherNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1Name, pt2Name, pt3Name))

bp <- ggplot(data = df,
             mapping = aes(x = Feature,
                           y = log2ObsExp,
                           fill = Sample)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("purple4",
                               "dodgerblue2",
                               "forestgreen"),
                    labels = c(pt1Name,
                               pt2Name,
                               pt3Name)) +
  geom_point(mapping = aes(Feature, log2alpha0.05),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey80", size = 10) +
  labs(y = expression("Log"[2]*"(observed/expected) peak overlap")) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  scale_x_discrete(position = "top") +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "top",
                             label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
                             nrow = 1,
                             byrow = TRUE)) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, angle = 90),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #legend.position = c(0.05, 0.30),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
        plot.title = element_text(size = 17, colour = "black", hjust = 0.5)) +
  ggtitle(paste0(plotTitle, " (", prettyNum(as.character(perms),
                                            big.mark = ",", trim = "T"),
                 " permutations)"))
ggsave(paste0(plotDir, "barplot_other_features_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_ASY1_CS_Rep1_ChIP_",
              pt1LibName, "_", pt2LibName, "_", pt3LibName, "_",
              region, "_peaks.pdf"),
       plot = bp,
       height = 8, width = 18)
save(bp,
     file = paste0(plotDir, "barplot_other_features_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_ASY1_CS_Rep1_ChIP_",
                   pt1LibName, "_", pt2LibName, "_", pt3LibName, "_",
                   region, "_peaks.RData"))


