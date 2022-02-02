#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 01.02.2022

# Examine relationships between differential fancm-wild type recombination rate (cM/Mb) and
# epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq

# Build GLM with differential fancm-wild type cM/Mb as the response variable
# and epigenetic and meiotic protein signals as predictor variables

# Usage:
# ./analyse_differential.R

options(stringsAsFactors = F)
library(fitdistrplus) # descdist, plotdist, fitdist included
library(glm2)
library(MASS) # glm.nb included; MASS is also loaded as a dependency of fitdistrplus
library(pscl) # zeroinfl included
library(vcd) # goodfit included
library(qualityTools) # qqPlot included
library(stats4) # mle (for estimating parameters by maximum likelihood) included
library(dplyr)
#library(segmentSeq)
#library(GenomicRanges)

#library(data.table)
#library(parallel)
#library(GenomicRanges)
#library(dplyr)
#library(plyr)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load cM/Mb data
dat <- read.table("AxC_mapped_marker_intervals_mean_ChIPseq_and_DNAmethyl_cMMb.tsv", header = T)
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
#[23] "wt_cM_inter"              "fancm_cM_inter"
#[25] "wt_cMMb_inter"            "fancm_cMMb_inter"
#[27] "diff_fancm_wt_cM_inter"   "diff_fancm_wt_cMMb_inter"
#[29] "divi_fancm_wt_cM_inter"   "divi_fancm_wt_cMMb_inter"
#[31] "l2fc_fancm_wt_cM_inter"   "l2fc_fancm_wt_cMMb_inter"

colnames(dat) <- c(colnames(dat)[1:9],
                   "ASY1", "DMC1", "H3K4me1", "H3K4me3", "H3K27ac",
                   "H3K27me3", "H3K36me3", "H3K9me2", "H3K27me1",
                   "CENH3", "mCG", "mCHG", "mCHH",
                   colnames(dat)[23:ncol(dat)])

pdf(paste0(plotDir, "hist_width_inter.pdf"))
hist(dat$width, breaks = 100)
dev.off()

pdf(paste0(plotDir, "hist_diff_fancm_wt_cM_inter.pdf"))
hist(dat$diff_fancm_wt_cM_inter, breaks = 100)
dev.off()

pdf(paste0(plotDir, "hist_diff_fancm_wt_cMMb_inter.pdf"))
hist(dat$diff_fancm_wt_cMMb_inter, breaks = 100)
dev.off()

pdf(paste0(plotDir, "hist_divi_fancm_wt_cM_inter.pdf"))
hist(dat$divi_fancm_wt_cM_inter, breaks = 100)
dev.off()

pdf(paste0(plotDir, "hist_divi_fancm_wt_cMMb_inter.pdf"))
hist(dat$divi_fancm_wt_cMMb_inter, breaks = 100)
dev.off()

pdf(paste0(plotDir, "hist_l2fc_fancm_wt_cM_inter.pdf"))
hist(dat$l2fc_fancm_wt_cM_inter, breaks = 100)
dev.off()

pdf(paste0(plotDir, "hist_l2fc_fancm_wt_cMMb_inter.pdf"))
hist(dat$l2fc_fancm_wt_cMMb_inter, breaks = 100)
dev.off()

# Inspect distribution of *_fancm_wt_cMMb_inter and *_fancm_wt_cM_inter:
# 1. by plotting the empirical density and the empirical cumulative distribution function (ECDF)
pdf(paste0(plotDir, "diff_fancm_wt_cMMb_inter_plotdist.pdf"))
fitdistrplus::plotdist(dat$diff_fancm_wt_cMMb_inter, histo = T, demp = T)
dev.off()

pdf(paste0(plotDir, "divi_fancm_wt_cMMb_inter_plotdist.pdf"))
fitdistrplus::plotdist(dat$divi_fancm_wt_cMMb_inter, histo = T, demp = T)
dev.off()

pdf(paste0(plotDir, "l2fc_fancm_wt_cMMb_inter_plotdist.pdf"))
fitdistrplus::plotdist(dat$l2fc_fancm_wt_cMMb_inter, histo = T, demp = T)
dev.off()

pdf(paste0(plotDir, "l2fc_fancm_wt_cM_inter_plotdist.pdf"))
fitdistrplus::plotdist(dat$l2fc_fancm_wt_cM_inter, histo = T, demp = T)
dev.off()

# Plot random normal distribution
set.seed(9382)
pdf(paste0(plotDir, "random_normal_n", length(dat$l2fc_fancm_wt_cMMb_inter), "_plotdist.pdf"))
fitdistrplus::plotdist(rnorm(n = length(dat$l2fc_fancm_wt_cMMb_inter),
                             mean = mean(dat$l2fc_fancm_wt_cMMb_inter),
                             sd = sd(dat$l2fc_fancm_wt_cMMb_inter)),
                       histo = T, demp = T)
dev.off()

# Plot random cauchy distribution
set.seed(9382)
pdf(paste0(plotDir, "random_cauchy_n", length(dat$l2fc_fancm_wt_cMMb_inter), "_plotdist.pdf"))
fitdistrplus::plotdist(rcauchy(n = length(dat$l2fc_fancm_wt_cMMb_inter)),
                       histo = T, demp = T)
dev.off()



# 2. by plotting a skewness-kurtosis (Cullen and Frey) graph
# using descdist from the fitdistrplus package, with 1000 bootstraps
# See https://stats.stackexchange.com/questions/58220/what-distribution-does-my-data-follow
# and https://stackoverflow.com/questions/31741742/how-to-identify-the-distribution-of-the-given-data-using-r
pdf(paste0(plotDir, "diff_fancm_wt_cMMb_inter_descdistPlot.pdf"))
fitdistrplus::descdist(dat$diff_fancm_wt_cMMb_inter, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  -1236.426   max:  1589.143
#median:  0
#mean:  4.930993
#estimated sd:  140.6925
#estimated skewness:  1.390143
#estimated kurtosis:  61.89126

pdf(paste0(plotDir, "divi_fancm_wt_cMMb_inter_descdistPlot.pdf"))
fitdistrplus::descdist(dat$divi_fancm_wt_cMMb_inter, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  0.0008081292   max:  1590.143
#median:  1
#mean:  16.2194
#estimated sd:  103.5066
#estimated skewness:  10.46197
#estimated kurtosis:  130.7225

pdf(paste0(plotDir, "l2fc_fancm_wt_cMMb_inter_descdistPlot.pdf"))
fitdistrplus::descdist(dat$l2fc_fancm_wt_cMMb_inter, discrete = F,
                       obs.col = "darkblue", boot = 1000, boot.col = "red")
dev.off()
#summary statistics
#------
#min:  -10.27313   max:  10.63494
#median:  0
#mean:  0.3096608
#estimated sd:  2.419111
#estimated skewness:  0.3494517
#estimated kurtosis:  8.749073



# Based on *_fancm_wt_cMMb_inter_descdistPlot.pdf, inspect fit to distributions
# Quantile-quantile (qq) plots for different distributiions
# cauchy
pdf(paste0(plotDir, "l2fc_fancm_wt_cMMb_inter_qqPlot_cauchy.pdf"))
qualityTools::qqPlot(x = dat$l2fc_fancm_wt_cMMb_inter, y = "cauchy")
dev.off()

pdf(paste0(plotDir, "l2fc_fancm_wt_cM_inter_qqPlot_cauchy.pdf"))
qualityTools::qqPlot(x = dat$l2fc_fancm_wt_cM_inter, y = "cauchy")
dev.off()

pdf(paste0(plotDir, "diff_fancm_wt_cMMb_inter_qqPlot_cauchy.pdf"))
qualityTools::qqPlot(x = dat$diff_fancm_wt_cMMb_inter, y = "cauchy")
dev.off()

pdf(paste0(plotDir, "diff_fancm_wt_cM_inter_qqPlot_cauchy.pdf"))
qualityTools::qqPlot(x = dat$diff_fancm_wt_cM_inter, y = "cauchy")
dev.off()

# normal
pdf(paste0(plotDir, "l2fc_fancm_wt_cMMb_inter_qqPlot_normal.pdf"))
qualityTools::qqPlot(x = dat$l2fc_fancm_wt_cMMb_inter, y = "normal")
dev.off()

pdf(paste0(plotDir, "l2fc_fancm_wt_cM_inter_qqPlot_normal.pdf"))
qualityTools::qqPlot(x = dat$l2fc_fancm_wt_cM_inter, y = "normal")
dev.off()

pdf(paste0(plotDir, "diff_fancm_wt_cMMb_inter_qqPlot_normal.pdf"))
qualityTools::qqPlot(x = dat$diff_fancm_wt_cMMb_inter, y = "normal")
dev.off()

pdf(paste0(plotDir, "diff_fancm_wt_cM_inter_qqPlot_normal.pdf"))
qualityTools::qqPlot(x = dat$diff_fancm_wt_cM_inter, y = "normal")
dev.off()



# Inspect fit of distribution using maximum likelihood estimation (MLE) within fitdistrplus1
# See https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
l2fc_fancm_wt_cMMb_inter_fitdist_cauchy <- fitdistrplus::fitdist(data = dat$l2fc_fancm_wt_cMMb_inter,
                                                                 distr = "cauchy",
                                                                 method = "mle")
l2fc_fancm_wt_cMMb_inter_fitdist_normal <- fitdistrplus::fitdist(data = dat$l2fc_fancm_wt_cMMb_inter,
                                                                 distr = "norm",
                                                                 method = "mle")

l2fc_fancm_wt_cM_inter_fitdist_cauchy <- fitdistrplus::fitdist(data = dat$l2fc_fancm_wt_cM_inter,
                                                               distr = "cauchy",
                                                               method = "mle")
l2fc_fancm_wt_cM_inter_fitdist_normal <- fitdistrplus::fitdist(data = dat$l2fc_fancm_wt_cM_inter,
                                                               distr = "norm",
                                                               method = "mle")


diff_fancm_wt_cMMb_inter_fitdist_cauchy <- fitdistrplus::fitdist(data = dat$diff_fancm_wt_cMMb_inter,
                                                                 distr = "cauchy",
                                                                 method = "mle")
diff_fancm_wt_cMMb_inter_fitdist_normal <- fitdistrplus::fitdist(data = dat$diff_fancm_wt_cMMb_inter,
                                                                 distr = "norm",
                                                                 method = "mle")

diff_fancm_wt_cM_inter_fitdist_cauchy <- fitdistrplus::fitdist(data = dat$diff_fancm_wt_cM_inter,
                                                               distr = "cauchy",
                                                               method = "mle")
diff_fancm_wt_cM_inter_fitdist_normal <- fitdistrplus::fitdist(data = dat$diff_fancm_wt_cM_inter,
                                                               distr = "norm",
                                                               method = "mle")


print(summary(l2fc_fancm_wt_cMMb_inter_fitdist_cauchy))
print(summary(l2fc_fancm_wt_cMMb_inter_fitdist_normal))
print(summary(l2fc_fancm_wt_cM_inter_fitdist_cauchy))
print(summary(l2fc_fancm_wt_cM_inter_fitdist_normal))

print(summary(diff_fancm_wt_cMMb_inter_fitdist_cauchy))
print(summary(diff_fancm_wt_cMMb_inter_fitdist_normal))
print(summary(diff_fancm_wt_cM_inter_fitdist_cauchy))
print(summary(diff_fancm_wt_cM_inter_fitdist_normal))


# Generate goodness-of-fit plots:
# See https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
# 1. "a density plot representing the density function of the fitted distribution
#     along with the histogram of the empirical distribution"
# 2. "a cumulative distribution function (CDF) plot of both the
#     empirical distribution and the fitted distribution"
# 3. "a Q-Q plot representing the empirical quantiles (y-axis)
#     against the theoretical quantiles (x-axis)"
# 4. "a P-P plot representing the empirical distribution function evaluated at each data point (y-axis)
#     against the fitted distribution function (x-axis)"

# Load modified denscomp and cdfcomp functions (allowing line widths to be changed)
load("fitdistrplus_denscomp_mod.RData")
load("fitdistrplus_cdfcomp_mod.RData")
pdf(paste0(plotDir, "l2fc_fancm_wt_cMMb_inter_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Cauchy", "Normal")
denscomp(list(l2fc_fancm_wt_cMMb_inter_fitdist_cauchy, l2fc_fancm_wt_cMMb_inter_fitdist_cauchy),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(l2fc_fancm_wt_cMMb_inter_fitdist_cauchy, l2fc_fancm_wt_cMMb_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(l2fc_fancm_wt_cMMb_inter_fitdist_cauchy, l2fc_fancm_wt_cMMb_inter_fitdist_cauchy),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(l2fc_fancm_wt_cMMb_inter_fitdist_cauchy, l2fc_fancm_wt_cMMb_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
dev.off()

pdf(paste0(plotDir, "l2fc_fancm_wt_cM_inter_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Cauchy", "Normal")
denscomp(list(l2fc_fancm_wt_cM_inter_fitdist_cauchy, l2fc_fancm_wt_cM_inter_fitdist_cauchy),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(l2fc_fancm_wt_cM_inter_fitdist_cauchy, l2fc_fancm_wt_cM_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(l2fc_fancm_wt_cM_inter_fitdist_cauchy, l2fc_fancm_wt_cM_inter_fitdist_cauchy),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(l2fc_fancm_wt_cM_inter_fitdist_cauchy, l2fc_fancm_wt_cM_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
dev.off()

pdf(paste0(plotDir, "diff_fancm_wt_cMMb_inter_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Cauchy", "Normal")
denscomp(list(diff_fancm_wt_cMMb_inter_fitdist_cauchy, diff_fancm_wt_cMMb_inter_fitdist_cauchy),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(diff_fancm_wt_cMMb_inter_fitdist_cauchy, diff_fancm_wt_cMMb_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(diff_fancm_wt_cMMb_inter_fitdist_cauchy, diff_fancm_wt_cMMb_inter_fitdist_cauchy),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(diff_fancm_wt_cMMb_inter_fitdist_cauchy, diff_fancm_wt_cMMb_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
dev.off()

pdf(paste0(plotDir, "diff_fancm_wt_cM_inter_GoFplots.pdf"))
par(mfrow = c(2, 2))
plot.legend <- c("Cauchy", "Normal")
denscomp(list(diff_fancm_wt_cM_inter_fitdist_cauchy, diff_fancm_wt_cM_inter_fitdist_cauchy),
         legendtext = plot.legend, fitlwd = 2)
qqcomp(list(diff_fancm_wt_cM_inter_fitdist_cauchy, diff_fancm_wt_cM_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
cdfcomp(list(diff_fancm_wt_cM_inter_fitdist_cauchy, diff_fancm_wt_cM_inter_fitdist_cauchy),
        legendtext = plot.legend, fitlwd = 2)
ppcomp(list(diff_fancm_wt_cM_inter_fitdist_cauchy, diff_fancm_wt_cM_inter_fitdist_cauchy),
       legendtext = plot.legend, fitpch = 20)
dev.off()


# First fit model using normal (Gaussian) distribution
glmNormal_l2fc_fancm_wt_cMMb_inter <- glm2(formula = l2fc_fancm_wt_cMMb_inter ~
                                                     (ASY1 + DMC1 + H3K4me1 + H3K4me3 + H3K27ac +
                                                      H3K27me3 + H3K36me3 + H3K9me2 + H3K27me1 +
                                                      CENH3 + mCG + mCHG + mCHH + width)^2,
                                           data = dat,
                                           family = gaussian(),
                                           control = glm.control(maxit = 100000))

summary(glmNormal_l2fc_fancm_wt_cMMb_inter)

glmNormal_l2fc_fancm_wt_cMMb_inter_stepAIC <- stepAIC(object = glmNormal_l2fc_fancm_wt_cMMb_inter, direction = "both")
print("stepAIC-selected model formula:")
print(formula(glmNormal_l2fc_fancm_wt_cMMb_inter_stepAIC))
#print(glmNormal_l2fc_fancm_wt_cMMb_inter_stepAIC$formula)
glmNormal_l2fc_fancm_wt_cMMb_inter_formula <- formula(glmNormal_l2fc_fancm_wt_cMMb_inter_stepAIC)
glmNormal_l2fc_fancm_wt_cMMb_inter_select <- glm2(formula = formula(glmNormal_l2fc_fancm_wt_cMMb_inter_stepAIC),
                        data = dat,
                        family = Gamma(link="inverse"),
                        control = glm.control(maxit = 100000))
stopifnot(identical(formula(glmNormal_l2fc_fancm_wt_cMMb_inter_select), formula(glmNormal_l2fc_fancm_wt_cMMb_inter_stepAIC)))
glmNormal_l2fc_fancm_wt_cMMb_inter_formula <- formula(glmNormal_l2fc_fancm_wt_cMMb_inter_select)
# Estimate the shape parameter and adjust GLM coefficient estimates and predictions
# See https://stats.stackexchange.com/questions/58497/using-r-for-glm-with-gamma-distribution
glmNormal_l2fc_fancm_wt_cMMb_inter_shape <- gamma.shape(glmNormal_l2fc_fancm_wt_cMMb_inter_select)
glmNormal_l2fc_fancm_wt_cMMb_inter_predict <- predict(glmNormal_l2fc_fancm_wt_cMMb_inter_select, type = "response",
                            se = T, dispersion = 1/glmNormal_l2fc_fancm_wt_cMMb_inter_shape$alpha)
glmNormal_l2fc_fancm_wt_cMMb_inter_summary <- summary(glmNormal_l2fc_fancm_wt_cMMb_inter_select, dispersion = 1/glmNormal_l2fc_fancm_wt_cMMb_inter_shape$alpha)
glmNormal_l2fc_fancm_wt_cMMb_inter_coeffs <- glmNormal_l2fc_fancm_wt_cMMb_inter_summary$coefficients



