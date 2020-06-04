#!/applications/R/R-3.5.0/bin/Rscript

# Plot bar graph of log2(observed:expected) peaks overlapping other features

# Usage:
# ./TEfams_vs_ASY1_euchromatin_heterochromatin_peaks_bargraphs_v030620.R "ASY1 peak overlap" Agenome Bgenome Dgenome Agenome 10000

library(ggplot2)
library(ggthemes)

#plotTitle <- "ASY1 peak overlap"
#pt1LibName <- "Agenome" 
#pt2LibName <- "Bgenome"
#pt3LibName <- "Dgenome"
#ptOrderLibName <- "Agenome"
## Number of permutations (randomisations) performed
#perms <- "10000"

args <- commandArgs(trailingOnly = T)
plotTitle <- args[1]
pt1LibName <- args[2]
pt2LibName <- args[3]
pt3LibName <- args[4]
ptOrderLibName <- args[5]
# Number of permutations (randomisations) performed
perms <- as.character(args[6])

ptOrderDir <- paste0("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                     "euchromatin/regioneR/noMinWidth_mergedOverlaps/")
euDir <- paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                "euchromatin/regioneR/noMinWidth_mergedOverlaps/")
heteroDir <- paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                    "heterochromatin/regioneR/noMinWidth_mergedOverlaps/")

plotDir <- "./bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

otherNamesPlot <- c(
                    "CACTA",
                    "Harbinger",
                    "hAT",
#4                    "Helitron",
                    "Mariner",
                    "Mutator",
                    "MITE",
                    "Unclass. class 2",
                    "Unclass. TIRs",
                    "Copia LTR",
                    "Gypsy LTR",
                    "LINE",
                    "SINE",
                    "Unclass. LTR",
                    "Unclass. repeats"
                   )

# Load permutation test results for peak set to be used for ordering
# of other features in bar graph
load(paste0(ptOrderDir, "permTest_", perms, "perms_DMC1_Rep1_ChIP_peaks_vs_TEfams_in_",
            ptOrderLibName, "_euchromatin.RData"))
ptOrder <- ptPeaksOtherPerChrom[-c(4)]
ptPeaksOtherPerChrom <- NULL

# Load permutation test results to be used for plotting
load(paste0(euDir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
            pt1LibName, "_euchromatin.RData"))
pt1eu <- ptPeaksOtherPerChrom[-c(4)]
ptPeaksOtherPerChrom <- NULL

load(paste0(euDir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
            pt2LibName, "_euchromatin.RData"))
pt2eu <- ptPeaksOtherPerChrom[-c(4)]
ptPeaksOtherPerChrom <- NULL

load(paste0(euDir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
            pt3LibName, "_euchromatin.RData"))
pt3eu <- ptPeaksOtherPerChrom[-c(4)]
ptPeaksOtherPerChrom <- NULL

load(paste0(heteroDir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
            pt1LibName, "_heterochromatin.RData"))
pt1hetero <- ptPeaksOtherPerChrom[-c(4)]
ptPeaksOtherPerChrom <- NULL

load(paste0(heteroDir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
            pt2LibName, "_heterochromatin.RData"))
pt2hetero <- ptPeaksOtherPerChrom[-c(4)]
ptPeaksOtherPerChrom <- NULL

load(paste0(heteroDir, "permTest_", perms, "perms_ASY1_CS_Rep1_ChIP_peaks_vs_TEfams_in_",
            pt3LibName, "_heterochromatin.RData"))
pt3hetero <- ptPeaksOtherPerChrom[-c(4)]
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
ptOrder_log2ObsExp_log2alpha0.05 <- unlist(ptOrder_log2ObsExp)-unlist(ptOrder_log2alpha0.05)

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

# pt1eu
pt1eu_Pval <- lapply(seq_along(pt1eu), function(x) {
  pt1eu[[x]]$numOverlaps$pval
})
pt1eu_Obs <- lapply(seq_along(pt1eu), function(x) {
  pt1eu[[x]]$numOverlaps$observed
})
pt1eu_Perm <- lapply(seq_along(pt1eu), function(x) {
  pt1eu[[x]]$numOverlaps$permuted
})
pt1eu_Exp <- lapply(seq_along(pt1eu), function(x) {
  mean(pt1eu[[x]]$numOverlaps$permuted)
})
pt1eu_log2ObsExp <- lapply(seq_along(pt1eu_Obs), function(x) {
  log2(pt1eu_Obs[[x]]/pt1eu_Exp[[x]])
})
pt1eu_Zscore <- lapply(seq_along(pt1eu), function(x) {
  pt1eu[[x]]$numOverlaps$zscore
})
pt1eu_AltHyp <- lapply(seq_along(pt1eu), function(x) {
  pt1eu[[x]]$numOverlaps$alternative
})
pt1eu_alpha0.05 <- lapply(seq_along(pt1eu_Perm), function(x) {
  if(pt1eu_AltHyp[[x]] == "greater") {
    quantile(pt1eu_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1eu_Perm[[x]], 0.05)[[1]]
  }
})
pt1eu_log2alpha0.05 <- lapply(seq_along(pt1eu_alpha0.05), function(x) {
  log2(pt1eu_alpha0.05[[x]]/pt1eu_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt1eu_log2ObsExp_sorted <- unlist(pt1eu_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt1eu_log2alpha0.05_sorted <- unlist(pt1eu_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt1eu_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt1eu_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                     decreasing = T,
                                                     index.return = T)$ix]

# pt2eu
pt2eu_Pval <- lapply(seq_along(pt2eu), function(x) {
  pt2eu[[x]]$numOverlaps$pval
})
pt2eu_Obs <- lapply(seq_along(pt2eu), function(x) {
  pt2eu[[x]]$numOverlaps$observed
})
pt2eu_Perm <- lapply(seq_along(pt2eu), function(x) {
  pt2eu[[x]]$numOverlaps$permuted
})
pt2eu_Exp <- lapply(seq_along(pt2eu), function(x) {
  mean(pt2eu[[x]]$numOverlaps$permuted)
})
pt2eu_log2ObsExp <- lapply(seq_along(pt2eu_Obs), function(x) {
  log2(pt2eu_Obs[[x]]/pt2eu_Exp[[x]])
})
pt2eu_Zscore <- lapply(seq_along(pt2eu), function(x) {
  pt2eu[[x]]$numOverlaps$zscore
})
pt2eu_AltHyp <- lapply(seq_along(pt2eu), function(x) {
  pt2eu[[x]]$numOverlaps$alternative
})
pt2eu_alpha0.05 <- lapply(seq_along(pt2eu_Perm), function(x) {
  if(pt2eu_AltHyp[[x]] == "greater") {
    quantile(pt2eu_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2eu_Perm[[x]], 0.05)[[1]]
  }
})
pt2eu_log2alpha0.05 <- lapply(seq_along(pt2eu_alpha0.05), function(x) {
  log2(pt2eu_alpha0.05[[x]]/pt2eu_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt2eu_log2ObsExp_sorted <- unlist(pt2eu_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt2eu_log2alpha0.05_sorted <- unlist(pt2eu_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt2eu_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt2eu_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# pt3eu
pt3eu_Pval <- lapply(seq_along(pt3eu), function(x) {
  pt3eu[[x]]$numOverlaps$pval
})
pt3eu_Obs <- lapply(seq_along(pt3eu), function(x) {
  pt3eu[[x]]$numOverlaps$observed
})
pt3eu_Perm <- lapply(seq_along(pt3eu), function(x) {
  pt3eu[[x]]$numOverlaps$permuted
})
pt3eu_Exp <- lapply(seq_along(pt3eu), function(x) {
  mean(pt3eu[[x]]$numOverlaps$permuted)
})
pt3eu_log2ObsExp <- lapply(seq_along(pt3eu_Obs), function(x) {
  log2(pt3eu_Obs[[x]]/pt3eu_Exp[[x]])
})
pt3eu_Zscore <- lapply(seq_along(pt3eu), function(x) {
  pt3eu[[x]]$numOverlaps$zscore
})
pt3eu_AltHyp <- lapply(seq_along(pt3eu), function(x) {
  pt3eu[[x]]$numOverlaps$alternative
})
pt3eu_alpha0.05 <- lapply(seq_along(pt3eu_Perm), function(x) {
  if(pt3eu_AltHyp[[x]] == "greater") {
    quantile(pt3eu_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt3eu_Perm[[x]], 0.05)[[1]]
  }
})
pt3eu_log2alpha0.05 <- lapply(seq_along(pt3eu_alpha0.05), function(x) {
  log2(pt3eu_alpha0.05[[x]]/pt3eu_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt3eu_log2ObsExp_sorted <- unlist(pt3eu_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt3eu_log2alpha0.05_sorted <- unlist(pt3eu_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt3eu_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt3eu_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# pt1hetero
pt1hetero_Pval <- lapply(seq_along(pt1hetero), function(x) {
  pt1hetero[[x]]$numOverlaps$pval
})
pt1hetero_Obs <- lapply(seq_along(pt1hetero), function(x) {
  pt1hetero[[x]]$numOverlaps$observed
})
pt1hetero_Perm <- lapply(seq_along(pt1hetero), function(x) {
  pt1hetero[[x]]$numOverlaps$permuted
})
pt1hetero_Exp <- lapply(seq_along(pt1hetero), function(x) {
  mean(pt1hetero[[x]]$numOverlaps$permuted)
})
pt1hetero_log2ObsExp <- lapply(seq_along(pt1hetero_Obs), function(x) {
  log2(pt1hetero_Obs[[x]]/pt1hetero_Exp[[x]])
})
pt1hetero_Zscore <- lapply(seq_along(pt1hetero), function(x) {
  pt1hetero[[x]]$numOverlaps$zscore
})
pt1hetero_AltHyp <- lapply(seq_along(pt1hetero), function(x) {
  pt1hetero[[x]]$numOverlaps$alternative
})
pt1hetero_alpha0.05 <- lapply(seq_along(pt1hetero_Perm), function(x) {
  if(pt1hetero_AltHyp[[x]] == "greater") {
    quantile(pt1hetero_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1hetero_Perm[[x]], 0.05)[[1]]
  }
})
pt1hetero_log2alpha0.05 <- lapply(seq_along(pt1hetero_alpha0.05), function(x) {
  log2(pt1hetero_alpha0.05[[x]]/pt1hetero_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt1hetero_log2ObsExp_sorted <- unlist(pt1hetero_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt1hetero_log2alpha0.05_sorted <- unlist(pt1hetero_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt1hetero_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt1hetero_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                     decreasing = T,
                                                     index.return = T)$ix]

# pt2hetero
pt2hetero_Pval <- lapply(seq_along(pt2hetero), function(x) {
  pt2hetero[[x]]$numOverlaps$pval
})
pt2hetero_Obs <- lapply(seq_along(pt2hetero), function(x) {
  pt2hetero[[x]]$numOverlaps$observed
})
pt2hetero_Perm <- lapply(seq_along(pt2hetero), function(x) {
  pt2hetero[[x]]$numOverlaps$permuted
})
pt2hetero_Exp <- lapply(seq_along(pt2hetero), function(x) {
  mean(pt2hetero[[x]]$numOverlaps$permuted)
})
pt2hetero_log2ObsExp <- lapply(seq_along(pt2hetero_Obs), function(x) {
  log2(pt2hetero_Obs[[x]]/pt2hetero_Exp[[x]])
})
pt2hetero_Zscore <- lapply(seq_along(pt2hetero), function(x) {
  pt2hetero[[x]]$numOverlaps$zscore
})
pt2hetero_AltHyp <- lapply(seq_along(pt2hetero), function(x) {
  pt2hetero[[x]]$numOverlaps$alternative
})
pt2hetero_alpha0.05 <- lapply(seq_along(pt2hetero_Perm), function(x) {
  if(pt2hetero_AltHyp[[x]] == "greater") {
    quantile(pt2hetero_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2hetero_Perm[[x]], 0.05)[[1]]
  }
})
pt2hetero_log2alpha0.05 <- lapply(seq_along(pt2hetero_alpha0.05), function(x) {
  log2(pt2hetero_alpha0.05[[x]]/pt2hetero_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt2hetero_log2ObsExp_sorted <- unlist(pt2hetero_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt2hetero_log2alpha0.05_sorted <- unlist(pt2hetero_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt2hetero_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt2hetero_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# pt3hetero
pt3hetero_Pval <- lapply(seq_along(pt3hetero), function(x) {
  pt3hetero[[x]]$numOverlaps$pval
})
pt3hetero_Obs <- lapply(seq_along(pt3hetero), function(x) {
  pt3hetero[[x]]$numOverlaps$observed
})
pt3hetero_Perm <- lapply(seq_along(pt3hetero), function(x) {
  pt3hetero[[x]]$numOverlaps$permuted
})
pt3hetero_Exp <- lapply(seq_along(pt3hetero), function(x) {
  mean(pt3hetero[[x]]$numOverlaps$permuted)
})
pt3hetero_log2ObsExp <- lapply(seq_along(pt3hetero_Obs), function(x) {
  log2(pt3hetero_Obs[[x]]/pt3hetero_Exp[[x]])
})
pt3hetero_Zscore <- lapply(seq_along(pt3hetero), function(x) {
  pt3hetero[[x]]$numOverlaps$zscore
})
pt3hetero_AltHyp <- lapply(seq_along(pt3hetero), function(x) {
  pt3hetero[[x]]$numOverlaps$alternative
})
pt3hetero_alpha0.05 <- lapply(seq_along(pt3hetero_Perm), function(x) {
  if(pt3hetero_AltHyp[[x]] == "greater") {
    quantile(pt3hetero_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt3hetero_Perm[[x]], 0.05)[[1]]
  }
})
pt3hetero_log2alpha0.05 <- lapply(seq_along(pt3hetero_alpha0.05), function(x) {
  log2(pt3hetero_alpha0.05[[x]]/pt3hetero_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt3hetero_log2ObsExp_sorted <- unlist(pt3hetero_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt3hetero_log2alpha0.05_sorted <- unlist(pt3hetero_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
#pt3hetero_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
#                                         decreasing = T,
#                                         index.return = T)$ix]
pt3hetero_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# Combine in data.frame
pt1Names <- paste0(substr(pt1LibName, 1, 1), "  ", c("R1 & R3", "R2a-R2b"))
pt2Names <- paste0(substr(pt2LibName, 1, 1), "  ", c("R1 & R3", "R2a-R2b"))
pt3Names <- paste0(substr(pt3LibName, 1, 1), "  ", c("R1 & R3", "R2a-R2b"))
pt1NamesPlot <- c("R1 & R3", "R2a-R2b")
pt2NamesPlot <- c("R1 & R3", "R2a-R2b")
pt3NamesPlot <- c("R1 & R3", "R2a-R2b")
df <- data.frame(Sample = rep(c(pt1Names, pt2Names, pt3Names),
                              each = length(ptOrder_log2ObsExp_sorted)),
                 plotSample = rep(c(pt1NamesPlot, pt2NamesPlot, pt3NamesPlot),
                                  each = length(ptOrder_log2ObsExp_sorted)),
                 Feature = rep(pt1eu_otherNamesPlot_sorted, 6),
                 log2ObsExp = c(pt1eu_log2ObsExp_sorted,
                                pt1hetero_log2ObsExp_sorted,
                                pt2eu_log2ObsExp_sorted,
                                pt2hetero_log2ObsExp_sorted,
                                pt3eu_log2ObsExp_sorted,
                                pt3hetero_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1eu_log2alpha0.05_sorted,
                                   pt1hetero_log2alpha0.05_sorted,
                                   pt2eu_log2alpha0.05_sorted,
                                   pt2hetero_log2alpha0.05_sorted,
                                   pt3eu_log2alpha0.05_sorted,
                                   pt3hetero_log2alpha0.05_sorted))

df$Feature <- factor(df$Feature,
                     levels = pt1eu_otherNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1Names, pt2Names, pt3Names))

bp <- ggplot(data = df,
             mapping = aes(x = Feature,
                           y = log2ObsExp,
                           fill = Sample)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("purple4", "purple1",
                               "dodgerblue4", "dodgerblue1",
                               "darkgreen", "limegreen"),
                    labels = c(pt1NamesPlot,
                               pt2NamesPlot,
                               pt3NamesPlot)) +
  geom_point(mapping = aes(Feature, log2alpha0.05),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey80", size = 14) +
  labs(y = expression("Log"[2]*"(observed/expected) overlap")) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  scale_x_discrete(position = "bottom") +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "right",
                             keywidth = 3, keyheight = 3,
                             label.theme = element_text(size = 44, hjust = 0, vjust = 0.5, angle = 0),
                             nrow = 3,
                             byrow = TRUE)) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 44, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 44, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 44, colour = "black", hjust = 0.5, vjust = 0.5, angle = 45),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.2, 0.2),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
        plot.title = element_text(size = 70, face = "bold", colour = "black", hjust = 0.5)) +
  ggtitle(paste0(plotTitle))#, " (", prettyNum(as.character(perms),
                 #                           big.mark = ",", trim = "T"),
                 #" permutations)"))
ggsave(paste0(plotDir, "barplot_TEfams_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_ASY1_CS_Rep1_ChIP_",
              pt1LibName, "_", pt2LibName, "_", pt3LibName,
              "_peaks_v030620.pdf"),
       plot = bp,
       height = 12, width = 36)
save(bp,
     file = paste0(plotDir, "barplot_TEfams_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_ASY1_CS_Rep1_ChIP_",
                   pt1LibName, "_", pt2LibName, "_", pt3LibName,
                   "_peaks_v030620.RData"))
