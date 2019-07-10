#!/applications/R/R-3.5.0/bin/Rscript

# Plot feature average SNP frequency profiles with 95% CIs

# Usage:
# /applications/R/R-3.5.0/bin/Rscript features_avgProfileRibbon.R ASY1_CS_peaks 'ASY1 peaks' 2000 2kb '2 kb' 20 20bp 'A' 'euchromatin' '60varieties' '60 varieties v. CS' 'deepskyblue4'

#featureName <- "ASY1_CS_peaks"
#featureNamePlot <- "ASY1 peaks"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#genomeName <- "A"
#region <- "euchromatin"
#profileNames <- unlist(strsplit("60varieties",
#                                split = ","))
#profileNamesPlot <- unlist(strsplit("60 varieties v. CS",
#                                    split = ","))
#colours <- unlist(strsplit("deepskyblue4",
#                           split = ","))

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

# Load feature SNP frequency matrix
featureMat <- as.matrix(read.table(paste0(matDir,
                                          featureName, "_in_",
                                          genomeName, "genome_", region,                             
                                          "_SNP_frequency_feature_smoothed_target_and_",
                                          flankName, "_flank_dataframe.txt"),
                                   header = T))

# Load ranLoc SNP frequency matrix
ranLocMat <- as.matrix(read.table(paste0(matDir,
                                         featureName, "_in_",
                                         genomeName, "genome_", region,                             
                                         "_SNP_frequency_ranLoc_smoothed_target_and_",
                                         flankName, "_flank_dataframe.txt"),
                                  header = T))

## feature
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFfeature <- data.frame(window = colnames(featureMat),
                            t(featureMat))

# Convert into tidy data.frame (long format)
tidyDFfeature <- gather(data  = wideDFfeature,
                        key   = feature,
                        value = coverage,
                        -window)

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
tidyDFfeature$window <- factor(tidyDFfeature$window,
                               levels = as.character(wideDFfeature$window))

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature <- data.frame(window = as.character(wideDFfeature$window),
                               n      = tapply(X     = tidyDFfeature$coverage,
                                               INDEX = tidyDFfeature$window,
                                               FUN   = length),
                               mean   = tapply(X     = tidyDFfeature$coverage,
                                               INDEX = tidyDFfeature$window,
                                               FUN   = mean,
                                               na.rm = T),
                               sd     = tapply(X     = tidyDFfeature$coverage,
                                               INDEX = tidyDFfeature$window,
                                               FUN   = sd,
                                               na.rm = T))

summaryDFfeature$window <- factor(summaryDFfeature$window,
                                  levels = as.character(wideDFfeature$window))
summaryDFfeature$winNo <- factor(1:dim(summaryDFfeature)[1])
summaryDFfeature$sem <- summaryDFfeature$sd/sqrt(summaryDFfeature$n-1)
summaryDFfeature$CI_lower <- summaryDFfeature$mean -
  qt(0.975, df = summaryDFfeature$n-1)*summaryDFfeature$sem
summaryDFfeature$CI_upper <- summaryDFfeature$mean +
  qt(0.975, df = summaryDFfeature$n-1)*summaryDFfeature$sem

summaryDFfeature$profileName <- profileNamesPlot

## ranLoc
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFranLoc <- data.frame(window = colnames(ranLocMat),
                           t(ranLocMat))

# Convert into tidy data.frame (long format)
tidyDFranLoc <- gather(data  = wideDFranLoc,
                       key   = ranLoc,
                       value = coverage,
                       -window)

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
tidyDFranLoc$window <- factor(tidyDFranLoc$window,
                              levels = as.character(wideDFranLoc$window))

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (ranLocs) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFranLoc <- data.frame(window = as.character(wideDFranLoc$window),
                              n      = tapply(X     = tidyDFranLoc$coverage,
                                              INDEX = tidyDFranLoc$window,
                                              FUN   = length),
                              mean   = tapply(X     = tidyDFranLoc$coverage,
                                              INDEX = tidyDFranLoc$window,
                                              FUN   = mean,
                                              na.rm = T),
                              sd     = tapply(X     = tidyDFranLoc$coverage,
                                              INDEX = tidyDFranLoc$window,
                                              FUN   = sd,
                                              na.rm = T))

summaryDFranLoc$window <- factor(summaryDFranLoc$window,
                                 levels = as.character(wideDFranLoc$window))
summaryDFranLoc$winNo <- factor(1:dim(summaryDFranLoc)[1])
summaryDFranLoc$sem <- summaryDFranLoc$sd/sqrt(summaryDFranLoc$n-1)
summaryDFranLoc$CI_lower <- summaryDFranLoc$mean -
  qt(0.975, df = summaryDFranLoc$n-1)*summaryDFranLoc$sem
summaryDFranLoc$CI_upper <- summaryDFranLoc$mean +
  qt(0.975, df = summaryDFranLoc$n-1)*summaryDFranLoc$sem

summaryDFranLoc$profileName <- profileNamesPlot

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

# Plot average coverage profiles with 95% CI ribbon
## feature
ggObjGA <- NULL
ggObj1 <- NULL
ggObj1 <- ggplot(data = summaryDFfeature,
                 mapping = aes(x = winNo,
                               y = mean,
                               group = profileName),
                ) +
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
                              (dim(summaryDFfeature)[1])-(downstream/binSize),
                              dim(summaryDFfeature)[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature)[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = expression("Windowed SNP frequency")) +
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
                              (dim(summaryDFranLoc)[1])-(downstream/binSize),
                              dim(summaryDFranLoc)[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranLoc)[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = "") +
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
