#!/applications/R/R-3.5.0/bin/Rscript

# Plot feature average SNP frequency profiles with 95% CIs

# Usage:
# /applications/R/R-3.5.0/bin/Rscript features_avgProfileRibbon_SNPsubsets2.R ASY1_CS_Rep1_ChIP_peaks 'ASY1 peaks' 2000 2kb '2 kb' 20 20bp 'A' 'euchromatin' 'upstream_gene_variant,downstream_gene_variant' 'Upstream,Downstream' 'blue,green' '0.02,0.96'

#'upstream_gene_variant,downstream_gene_variant'
#'HIGH,MODERATE,MODIFIER,LOW'
#'missense_variant,synonymous_variant,transversion,transition'
#all
#upstream_gene_variant
#downstream_gene_variant
#HIGH
#MODERATE
#MODIFIER
#LOW
#missense_variant
#synonymous_variant
#transversion
#transition

#featureName <- "ASY1_CS_Rep1_ChIP_peaks"
#featureNamePlot <- "ASY1 peaks"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#genomeName <- "A"
#region <- "euchromatin"
#profileNames <- unlist(strsplit("HIGH,MODERATE,MODIFIER,LOW",
#                                split = ","))
#profileNamesPlot <- unlist(strsplit("High,Moderate,Modifier,Low",
#                                    split = ","))
#colours <- unlist(strsplit("red,orange,green,darkgreen",
#                           split = ","))
### Custom annotation legends using textGrob()
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## middle left upper
#legendPos <- as.numeric(unlist(strsplit("0.02,0.70",
#                                        split = ",")))
## middle left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.52",
#                                        split = ",")))
## middle left lower
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
## middle right
#legendPos <- as.numeric(unlist(strsplit("0.70,0.54",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
#                                        split = ",")))

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
featureNamePlot <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
genomeName <- args[8]
region <- args[9]
profileNames <- unlist(strsplit(args[10],
                                split = ","))
profileNamesPlot <- unlist(strsplit(args[11],
                                    split = ","))
colours <- unlist(strsplit(args[12],
                           split = ","))
legendPos <- as.numeric(unlist(strsplit(args[13],
                                        split = ",")))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

matDir <- "matrices/"
plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load feature SNP frequency matrix for each SNP subset
featureMats <- mclapply(seq_along(profileNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              featureName, "_in_",
                              genomeName, "genome_", region, "_",
                              profileNames[x],
                              "_SNP_frequency_feature_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(profileNames))

# Load ranLoc SNP frequency matrix for each SNP subset
ranLocMats <- mclapply(seq_along(profileNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              featureName, "_in_",
                              genomeName, "genome_", region, "_",
                              profileNames[x],
                              "_SNP_frequency_ranLoc_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(profileNames))

## feature
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFfeature_list <- mclapply(seq_along(featureMats), function(x) {
  data.frame(window = colnames(featureMats[[x]]),
             t(featureMats[[x]]))
}, mc.cores = length(featureMats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list  <- mclapply(seq_along(wideDFfeature_list), function(x) {
  gather(data  = wideDFfeature_list[[x]],
         key   = feature,
         value = coverage,
         -window)
}, mc.cores = length(wideDFfeature_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list)) {
  tidyDFfeature_list[[x]]$window <- factor(tidyDFfeature_list[[x]]$window,
                                           levels = as.character(wideDFfeature_list[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list  <- mclapply(seq_along(tidyDFfeature_list), function(x) {
  data.frame(window = as.character(wideDFfeature_list[[x]]$window),
             n      = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = mean,
                             na.rm = T),
             sd     = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = sd,
                             na.rm = T))
}, mc.cores = length(tidyDFfeature_list))

for(x in seq_along(summaryDFfeature_list)) {
  summaryDFfeature_list[[x]]$window <- factor(summaryDFfeature_list[[x]]$window,
                                              levels = as.character(wideDFfeature_list[[x]]$window))
  summaryDFfeature_list[[x]]$winNo <- factor(1:dim(summaryDFfeature_list[[x]])[1])
  summaryDFfeature_list[[x]]$sem <- summaryDFfeature_list[[x]]$sd/sqrt(summaryDFfeature_list[[x]]$n-1)
  summaryDFfeature_list[[x]]$CI_lower <- summaryDFfeature_list[[x]]$mean -
    qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
  summaryDFfeature_list[[x]]$CI_upper <- summaryDFfeature_list[[x]]$mean +
    qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
}

names(summaryDFfeature_list) <- profileNamesPlot

# Convert list summaryDFfeature_list into a single data.frame for plotting
summaryDFfeature <- bind_rows(summaryDFfeature_list, .id = "profileName")
summaryDFfeature$profileName <- factor(summaryDFfeature$profileName,
                                       levels = names(summaryDFfeature_list))

## ranLoc
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFranLoc_list <- mclapply(seq_along(ranLocMats), function(x) {
  data.frame(window = colnames(ranLocMats[[x]]),
             t(ranLocMats[[x]]))
}, mc.cores = length(ranLocMats))

# Convert into tidy data.frame (long format)
tidyDFranLoc_list  <- mclapply(seq_along(wideDFranLoc_list), function(x) {
  gather(data  = wideDFranLoc_list[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFranLoc_list)) {
  tidyDFranLoc_list[[x]]$window <- factor(tidyDFranLoc_list[[x]]$window,
                                          levels = as.character(wideDFranLoc_list[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (ranLocs) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFranLoc_list  <- mclapply(seq_along(tidyDFranLoc_list), function(x) {
  data.frame(window = as.character(wideDFranLoc_list[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = mean,
                             na.rm = T),
             sd     = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = sd,
                             na.rm = T))
}, mc.cores = length(tidyDFranLoc_list))

for(x in seq_along(summaryDFranLoc_list)) {
  summaryDFranLoc_list[[x]]$window <- factor(summaryDFranLoc_list[[x]]$window,
                                             levels = as.character(wideDFranLoc_list[[x]]$window))
  summaryDFranLoc_list[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list[[x]])[1])
  summaryDFranLoc_list[[x]]$sem <- summaryDFranLoc_list[[x]]$sd/sqrt(summaryDFranLoc_list[[x]]$n-1)
  summaryDFranLoc_list[[x]]$CI_lower <- summaryDFranLoc_list[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
  summaryDFranLoc_list[[x]]$CI_upper <- summaryDFranLoc_list[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
}

names(summaryDFranLoc_list) <- profileNamesPlot

# Convert list summaryDFranLoc_list into a single data.frame for plotting
summaryDFranLoc <- bind_rows(summaryDFranLoc_list, .id = "profileName")
summaryDFranLoc$profileName <- factor(summaryDFranLoc$profileName,
                                      levels = names(summaryDFranLoc_list))


if(featureName == "Genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Define y-axis limits
#ymin <- min(c(summaryDFfeature$mean-summaryDFfeature$sem,
#              summaryDFranLoc$mean-summaryDFranLoc$sem))
#ymax <- max(c(summaryDFfeature$mean+summaryDFfeature$sem,
#              summaryDFranLoc$mean+summaryDFranLoc$sem))
ymin <- min(c(summaryDFfeature$CI_lower,
              summaryDFranLoc$CI_lower))
ymax <- max(c(summaryDFfeature$CI_upper,
              summaryDFranLoc$CI_upper))

# Function for formatting y-axis labels
# with a given number of decimals
fmt_decimals <- function(decimals) {
  function(x) format(x, nsmall = decimals, scientific = FALSE)
}

# Define legend labels
legendLabs <- lapply(seq_along(profileNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(profileNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = colours[x], fontsize = 18)))
})

# Plot average coverage profiles with 95% CI ribbon
## feature
ggObjGA <- NULL
ggObj1 <- NULL
ggObj1 <- ggplot(data = summaryDFfeature,
                 mapping = aes(x = winNo,
                               y = mean,
                               group = profileName)) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = profileName),
            size = 1) +
  scale_colour_manual(values = colours) +
  geom_ribbon(data = summaryDFfeature,
              #mapping = aes(ymin = mean-sem,
              #              ymax = mean+sem,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = profileName),
              alpha = 0.4) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(ymin, ymax),
                     labels = function(x) sprintf("%7.5f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_list[[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = expression("Windowed SNP frequency")) +
  annotation_custom(legendLabs[[1]]) +
  annotation_custom(legendLabs[[2]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        #legend.text = element_text(size = 10),
        #legend.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(colour = "transparent",
        #                          fill = "transparent"),
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(.(featureNamePlot) ~ "in" ~ .(genomeName)*"-genome" ~ .(region) ~
                 "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
## ranLoc
ggObj2 <- NULL
ggObj2 <- ggplot(data = summaryDFranLoc,
                 mapping = aes(x = winNo,
                               y = mean,
                               group = profileName),
                ) +
  geom_line(data = summaryDFranLoc,
            mapping = aes(colour = profileName),
            size = 1) +
  scale_colour_manual(values = colours) +
  geom_ribbon(data = summaryDFranLoc,
              #mapping = aes(ymin = mean-sem,
              #              ymax = mean+sem,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = profileName),
              alpha = 0.4) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(ymin, ymax),
                     labels = function(x) sprintf("%7.5f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranLoc_list[[1]])[1])-(downstream/binSize),
                              dim(summaryDFranLoc_list[[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranLoc_list[[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = "") +
  annotation_custom(legendLabs[[1]]) +
  annotation_custom(legendLabs[[2]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        #legend.text = element_text(size = 10),
        #legend.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(colour = "transparent",
        #                          fill = "transparent"),
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote("Random loci in" ~ .(genomeName)*"-genome" ~ .(region) ~
                 "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranLoc$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
ggObjGA <- grid.arrange(ggObj1, ggObj2, nrow = 1, ncol = 2)
ggsave(paste0(plotDir, "SNPfreq_",
              paste(profileNames, collapse = "_"),
              "_around_", featureName, "_in_",
              genomeName, "genome_", region, ".pdf"),
       plot = ggObjGA,
       height = 7, width = 16)
