#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) peaks overlapping other features

# Usage:
# /applications/R/R-3.3.2/bin/Rscript TE_family_vs_H3Kmod_genome_wide_peaks_chr.R "Histone H3 lysine modification peaks" 10000 chr1A

library(ggplot2)
library(ggthemes)

dataName <- "Histone H3 lysine modification peaks"
perms <- 10000
chrName <- "chr1A"

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[2])
chrName <- args[3]

plotDir <- paste0("./", chrName, "/plots/")

ptDirs <- c(
            paste0("/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/", chrName, "/"),
            paste0("/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/", chrName, "/"),
#            paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K4me3/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/", chrName, "/"),
            paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K9ac/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/", chrName, "/"),
            paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/", chrName, "/"),
            paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K36me3/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/", chrName, "/")
           )

ptLibNames <- c(
                "H3K9me2_Rep1_ChIP",
                "H3K4me3_Rep1_ChIP",
#                "H3K4me3_ChIP_SRR6350668",
                "H3K9ac_ChIP_SRR6350667",
                "H3K27me3_ChIP_SRR6350666",
                "H3K36me3_ChIP_SRR6350670"
               )
ptLibNamesPlot <- c(
                    "H3K9me2",
                    "H3K4me3",
#                    "H3K4me3 (IWGSC)",
                    "H3K9ac",
                    "H3K27me3",
                    "H3K36me3"
                   )

famNames <- c(
              "CACTA_DTC",
              "Harbinger_DTH",
              "hAT_DTA",
              "Helitrons_DHH",
              "Mariner_DTT",
              "Mutator_DTM",
              "Unclassified_class_2_DXX",
              "Unclassified_with_TIRs_DTX",
              "Copia_LTR_RLC",
              "Gypsy_LTR_RLG",
              "LINE_RIX",
              "SINE_SIX",
              "Unclassified_LTR_RLX",
              "Unclassified_repeats_XXX"
             )
famNamesPlot <- c(
                  "CACTA",
                  "Harbinger",
                  "hAT",
                  "Helitrons",
                  "Mariner",
                  "Mutator",
                  "Unclass. class 2",
                  "Unclass. with TIRs",
                  "Copia LTR",
                  "Gypsy LTR",
                  "LINE",
                  "SINE",
                  "Unclass. LTR",
                  "Unclass. repeats"
                 )

ptList <- lapply(seq_along(ptLibNames), function(x) {
  load(paste0(ptDirs[x], "permTest_", ptLibNames[x],
              "_rangerPeaks_vs_TEsDNA_", chrName, ".RData"))
  load(paste0(ptDirs[x], "permTest_", ptLibNames[x],
              "_rangerPeaks_vs_TEsRNA_", chrName, ".RData"))
  c(ptPeaksTEsDNAPerChrom, ptPeaksTEsRNAPerChrom)
})
ptPeaksTEsDNAPerChrom <- NULL
ptPeaksTEsRNAPerChrom <- NULL

#assign(paste0(ptLibNames[2], "_TEsDNA"), ptPeaksTEsDNAPerChrom)

ptList_Pval <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList[[x]]), function(y) {
    ptList[[x]][[y]]$numOverlaps$pval
  })
})
ptList_Obs <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList[[x]]), function(y) {
    ptList[[x]][[y]]$numOverlaps$observed
  })
})
ptList_Perm <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList[[x]]), function(y) {
    ptList[[x]][[y]]$numOverlaps$permuted
  })
})
ptList_Exp <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList[[x]]), function(y) {
    mean(ptList[[x]][[y]]$numOverlaps$permuted)
  })
})
ptList_log2ObsExp <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList[[x]]), function(y) {
    log2((ptList_Obs[[x]][[y]]+1)/(ptList_Exp[[x]][[y]]+1))
  })
})
ptList_Zscore <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList[[x]]), function(y) {
    ptList[[x]][[y]]$numOverlaps$zscore
  })
})
ptList_AltHyp <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList[[x]]), function(y) {
    ptList[[x]][[y]]$numOverlaps$alternative
  })
})
ptList_alpha0.05 <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList_Perm[[x]]), function(y) {
    if(ptList_AltHyp[[x]][[y]] == "greater") {
      quantile(ptList_Perm[[x]][[y]], 0.95)[[1]]
    } else {
      quantile(ptList_Perm[[x]][[y]], 0.05)[[1]]
    }
  })
})
ptList_log2alpha0.05 <- lapply(seq_along(ptList), function(x) {
  lapply(seq_along(ptList_alpha0.05[[x]]), function(y) {
    log2((ptList_alpha0.05[[x]][[y]]+1)/(ptList_Exp[[x]][[y]]+1))
  })
})

ptList_log2ObsExp_sorted <- lapply(seq_along(ptList), function(x) {
  unlist(ptList_log2ObsExp[[x]][sort.int(unlist(ptList_log2ObsExp[[1]]),
                                         decreasing = T,
                                         index.return = T)$ix])
})
ptList_log2alpha0.05_sorted <- lapply(seq_along(ptList), function(x) {
  unlist(ptList_log2alpha0.05[[x]][sort.int(unlist(ptList_log2ObsExp[[1]]),
                                            decreasing = T,
                                            index.return = T)$ix])
})
famNames_sorted <- famNames[sort.int(unlist(ptList_log2ObsExp[[1]]),
                                     decreasing = T,
                                     index.return = T)$ix]
famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(ptList_log2ObsExp[[1]]),
                                             decreasing = T,
                                             index.return = T)$ix]

df <- data.frame(Sample = rep(ptLibNames,
                              each = length(ptList_log2ObsExp_sorted[[1]])),
                 Transposon_family = rep(famNamesPlot_sorted, length(ptLibNames)),
                 log2ObsExp = c(ptList_log2ObsExp_sorted[[1]],
                                ptList_log2ObsExp_sorted[[2]],
                                ptList_log2ObsExp_sorted[[3]],
                                ptList_log2ObsExp_sorted[[4]],
                                ptList_log2ObsExp_sorted[[5]]),
                 log2alpha0.05 = c(ptList_log2alpha0.05_sorted[[1]],
                                   ptList_log2alpha0.05_sorted[[2]],
                                   ptList_log2alpha0.05_sorted[[3]],
                                   ptList_log2alpha0.05_sorted[[4]],
                                   ptList_log2alpha0.05_sorted[[5]]))

df$Transposon_family <- factor(df$Transposon_family,
                               levels = famNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = ptLibNames)

bp <- ggplot(data = df,
             mapping = aes(x = Transposon_family,
                           y = log2ObsExp,
                           fill = Sample)) +
      theme_bw(base_size = 20) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Sample",
                        values = c("magenta3",
                                   "forestgreen",
#                                   "green2",
                                   "dodgerblue",
                                   "navy",
                                   "darkorange2"),
                        labels = ptLibNamesPlot) +
      geom_point(mapping = aes(Transposon_family, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "grey70", size = 7) +
      labs(x = "Transposon superfamily",
           y = expression("Log"[2]*"(observed:expected) peak overlap")) +
      theme(axis.line.y = element_line(size = 1, colour = "black"),
            axis.ticks.y = element_line(size = 1, colour = "black"),
            axis.text.y = element_text(colour = "black", size = 18),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 18),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 20)) +
      ggtitle(paste0(dataName, " (", chrName, "; ", as.character(perms), " permutations)"))
ggsave(paste0(plotDir, "barplot_TE_family_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_H3Kmod_peaks.pdf"),
       plot = bp,
       height = 14, width = 20)
save(bp,
     file = paste0(plotDir, "barplot_TE_family_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_H3Kmod_peaks.RData"))
