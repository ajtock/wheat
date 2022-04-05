#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 04.04.2022

# Examine relationships between differential fancm-wild type recombination rate (cM/Mb) and
# epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq

# Build GLM with differential fancm-wild type cM/Mb as the response variable
# and epigenetic and meiotic protein signals as predictor variables

# Usage:
# ./analyse_differential_cMMb_multicollinearity.R

options(stringsAsFactors = F)
library(glmnet)
library(caret)
library(dplyr)
library(glm2)
library(MASS) # stepAIC and glm.nb included; MASS is also loaded as a dependency of fitdistrplus
library(vcd) # goodfit included
library(car) # vif included
#library(fitdistrplus) # descdist, plotdist, fitdist included
#library(pscl) # zeroinfl included
#library(qualityTools) # qqPlot included
#library(stats4) # mle (for estimating parameters by maximum likelihood) included
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

coefDir <- "coefficients/"
plotDir <- "plots/"
system(paste0("[ -d ", coefDir, " ] || mkdir ", coefDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load response and predictors data
dat <- read.table("../../AxC_mapped_marker_intervals_mean_ChIPseq_DNAmethyl_cMMb_distanceToCEN.tsv", header = T)
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
#[25] "distToCEN"                "wt_cM_inter"             
#[27] "fancm_cM_inter"           "wt_cMMb_inter"           
#[29] "fancm_cMMb_inter"         "diff_fancm_wt_cM_inter"  
#[31] "diff_fancm_wt_cMMb_inter" "divi_fancm_wt_cM_inter"  
#[33] "divi_fancm_wt_cMMb_inter" "l2fc_fancm_wt_cM_inter"  
#[35] "l2fc_fancm_wt_cMMb_inter"

colnames(dat) <- c(colnames(dat)[1:9],
                   "ASY1", "DMC1", "H3K4me1", "H3K4me3", "H3K27ac",
                   "H3K27me3", "H3K36me3", "H3K9me2", "H3K27me1",
                   "CENH3", "mCG", "mCHG", "mCHH", "cMMb_1Mb", "cMMb_10Mb", "distToCEN",
                   colnames(dat)[26:ncol(dat)])

# Remove "cMMb_1Mb" and "distToCEN"
dat <- dat[,which(colnames(dat) %in% c("chr", "start", "end", "diff_fancm_wt_cMMb_inter",
                                       "ASY1", "DMC1", "H3K4me1", "H3K4me3", "H3K27ac",
                                       "H3K27me3", "H3K36me3", "H3K9me2", "H3K27me1",
#                                       "CENH3", "mCG", "mCHG", "mCHH"))]
                                       "CENH3", "mCG", "mCHG", "mCHH", "cMMb_10Mb"))]

colnames(dat)[which(colnames(dat) == "diff_fancm_wt_cMMb_inter")] <- "Diff_cMMb"
colnames(dat)[which(colnames(dat) == "cMMb_10Mb")] <- "IWGSC_cMMb"


# Identify predictor variables that show multicollinearity

# First centre and standardize the predictor variables to reduce
# structural multicollinearity (i.e., resulting from interaction terms)
# https://statisticsbyjim.com/regression/multicollinearity-in-regression-analysis/
# See https://advstats.psychstat.org/book/mregression/standardization.php
# https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
x_cols <- colnames(dat)[which(!(colnames(dat) %in% c("Diff_cMMb", "chr", "start", "end")))]
pre_proc_val <- caret::preProcess(dat[,x_cols], method = c("center", "scale"))
dat_scaled <- dat
dat_scaled[,x_cols] <- stats::predict(pre_proc_val, newdata = dat_scaled[,x_cols])
summary(dat_scaled)

## glmNormal_Diff_cMMb
# First fit model using normal (Gaussian) distribution
glmNormal_Diff_cMMb <- glm2(formula = Diff_cMMb ~
                                      (IWGSC_cMMb + ASY1 + DMC1 + H3K4me1 + H3K4me3 + H3K27ac +
                                       H3K27me3 + H3K36me3 + H3K9me2 + H3K27me1 +
                                       CENH3 + mCG + mCHG + mCHH),
                            data = dat_scaled,
                            family = gaussian(),
                            control = glm.control(maxit = 100000))
summary(glmNormal_Diff_cMMb)

# Use car::vif() to calculate generalized variance-inflation factors (GVIFs)
# for predictor variables
# https://statisticsbyjim.com/regression/multicollinearity-in-regression-analysis/
# https://www.statology.org/variance-inflation-factor-r/
# "VIFs greater than 5 represent critical levels of multicollinearity
# where the coefficients are poorly estimated, and the p-values are questionable."
# "A value greater than 5 indicates potentially severe correlation between a given
# predictor variable and other predictor variables in the model. In this case,
# the coefficient estimates and p-values in the regression output are likely unreliable."
glmNormal_Diff_cMMb_vif <- car::vif(mod = glmNormal_Diff_cMMb)
print(glmNormal_Diff_cMMb_vif)
#IWGSC_cMMb       ASY1       DMC1    H3K4me1    H3K4me3    H3K27ac   H3K27me3                                                                  
#  1.404341   4.659480   4.232534   4.052755  10.794569   2.382827   3.982482                                                                  
#  H3K36me3    H3K9me2   H3K27me1      CENH3        mCG       mCHG       mCHH                                                                  
#  7.948335   9.796587   6.544193   1.347513   4.521296   9.716677   1.458983
# First remove predictors with VIF > 7
glmNormal_Diff_cMMb_vif_gt7 <- names(which(glmNormal_Diff_cMMb_vif > 7))
x_retained <- names(dat_scaled)[ which( !( names(dat_scaled) %in%
                                           c(glmNormal_Diff_cMMb_vif_gt7, "Diff_cMMb", "chr", "start", "end") ) ) ]
print(x_retained)
#  [1] "ASY1"       "DMC1"       "H3K4me1"    "H3K27ac"    "H3K27me3"  
# [6] "H3K27me1"   "CENH3"      "mCG"        "mCHH"       "IWGSC_cMMb"

# Fit model again, this time using predictors with VIF <= 7
# Also remove ASY1 due to strong correlation with DMC1 and higher VIF
glmNormal_Diff_cMMb <- glm2(formula = Diff_cMMb ~
                                      (IWGSC_cMMb + DMC1 + H3K4me1 + H3K27ac +
                                       H3K27me3 + H3K27me1 + CENH3 + mCG + mCHH)^2,
                            data = dat_scaled,
                            family = gaussian(),
                            control = glm.control(maxit = 100000))
summary(glmNormal_Diff_cMMb)
#glmNormal_Diff_cMMb_vif <- car::vif(mod = glmNormal_Diff_cMMb)
#print(glmNormal_Diff_cMMb_vif)

glmNormal_Diff_cMMb_stepAIC <- MASS::stepAIC(object = glmNormal_Diff_cMMb, direction = "both")
print("stepAIC-selected model formula:")
print(formula(glmNormal_Diff_cMMb_stepAIC))
#print(glmNormal_Diff_cMMb_stepAIC$formula)
glmNormal_Diff_cMMb_formula <- formula(glmNormal_Diff_cMMb_stepAIC)
glmNormal_Diff_cMMb_select <- glm2(formula = formula(glmNormal_Diff_cMMb_stepAIC),
                                   data = dat_scaled,
                                   family = gaussian(),
                                   control = glm.control(maxit = 100000))
stopifnot(identical(formula(glmNormal_Diff_cMMb_select), formula(glmNormal_Diff_cMMb_stepAIC)))
glmNormal_Diff_cMMb_summary <- summary(glmNormal_Diff_cMMb_select)
glmNormal_Diff_cMMb_coeffs <- glmNormal_Diff_cMMb_summary$coefficients

glmNormal_Diff_cMMb_vif <- car::vif(mod = glmNormal_Diff_cMMb_select)
print(glmNormal_Diff_cMMb_vif)
#         IWGSC_cMMb                DMC1             H3K4me1             H3K27ac 
#           2.130386            3.384301            5.258543            9.103506 
#           H3K27me3            H3K27me1               CENH3                 mCG 
#           9.335585            6.085337            4.022470            6.984738 
#               mCHH     IWGSC_cMMb:DMC1  IWGSC_cMMb:H3K4me1  IWGSC_cMMb:H3K27ac 
#           5.117250            1.414272            2.895916            8.827602 
#IWGSC_cMMb:H3K27me3 IWGSC_cMMb:H3K27me1    IWGSC_cMMb:CENH3      IWGSC_cMMb:mCG 
#           2.959334            3.141913            7.550011            3.467106 
#       DMC1:H3K4me1        DMC1:H3K27ac          DMC1:CENH3           DMC1:mCHH 
#           3.133825            6.339808            3.793381            2.357738 
#    H3K4me1:H3K27ac    H3K4me1:H3K27me3       H3K4me1:CENH3         H3K4me1:mCG 
#           5.930169            4.824775            3.632953            4.581755 
#       H3K4me1:mCHH    H3K27ac:H3K27me3    H3K27ac:H3K27me1       H3K27ac:CENH3 
#           3.121945            2.306274           10.371693            5.565287 
#       H3K27ac:mCHH   H3K27me3:H3K27me1      H3K27me3:CENH3        H3K27me3:mCG 
#           2.991389            4.731663            5.704902            3.372501 
#      H3K27me3:mCHH      H3K27me1:CENH3        H3K27me1:mCG          CENH3:mCHH 
#           2.125289            5.835996            3.845523            2.474612 
#           mCG:mCHH 
#           3.705559

# Remove interaction terms with VIF > or ~ 5 (H3K4me1:H3K27me3 and H3K27me3:H3K27me1 are ~ 5)
glmNormal_Diff_cMMb <- glm2(formula = Diff_cMMb ~
                                      (IWGSC_cMMb + DMC1 + H3K4me1 + H3K27ac + H3K27me3 +
                                       H3K27me1 + CENH3 + mCG + mCHH +
                                       IWGSC_cMMb:DMC1 + IWGSC_cMMb:H3K4me1 + IWGSC_cMMb:H3K27me3 +
                                       IWGSC_cMMb:H3K27me1 + IWGSC_cMMb:mCG +
                                       DMC1:H3K4me1 + DMC1:CENH3 + DMC1:mCHH +
                                       H3K4me1:CENH3 + H3K4me1:mCG + H3K4me1:mCHH +
                                       H3K27ac:H3K27me3 + H3K27ac:mCHH +
                                       H3K27me3:mCG + H3K27me3:mCHH +
                                       H3K27me1:mCG + CENH3:mCHH + mCG:mCHH),
                            data = dat_scaled,
                            family = gaussian(),
                            control = glm.control(maxit = 100000))
summary(glmNormal_Diff_cMMb)
glmNormal_Diff_cMMb_vif <- car::vif(mod = glmNormal_Diff_cMMb)
print(glmNormal_Diff_cMMb_vif)
#         IWGSC_cMMb                DMC1             H3K4me1             H3K27ac 
#           1.800531            3.068502            3.449827            2.561653 
#           H3K27me3            H3K27me1               CENH3                 mCG 
#           3.950694            4.080383            1.708811            5.117890 
#               mCHH     IWGSC_cMMb:DMC1  IWGSC_cMMb:H3K4me1 IWGSC_cMMb:H3K27me3 
#           4.146239            1.290076            1.467827            2.393176 
#IWGSC_cMMb:H3K27me1      IWGSC_cMMb:mCG        DMC1:H3K4me1          DMC1:CENH3 
#           2.221078            1.927802            2.645060            2.398202 
#          DMC1:mCHH       H3K4me1:CENH3         H3K4me1:mCG        H3K4me1:mCHH 
#           2.243012            2.571169            3.262969            2.916701 
#   H3K27ac:H3K27me3        H3K27ac:mCHH        H3K27me3:mCG       H3K27me3:mCHH 
#           1.820626            2.168158            2.718596            1.774556 
#       H3K27me1:mCG          CENH3:mCHH            mCG:mCHH 
#           2.942261            2.059523            3.148928 
                                       
glmNormal_Diff_cMMb_stepAIC <- MASS::stepAIC(object = glmNormal_Diff_cMMb, direction = "both")
print("stepAIC-selected model formula:")
print(formula(glmNormal_Diff_cMMb_stepAIC))
#print(glmNormal_Diff_cMMb_stepAIC$formula)
glmNormal_Diff_cMMb_formula <- formula(glmNormal_Diff_cMMb_stepAIC)
glmNormal_Diff_cMMb_select <- glm2(formula = formula(glmNormal_Diff_cMMb_stepAIC),
                                   data = dat_scaled,
                                   family = gaussian(),
                                   control = glm.control(maxit = 100000))
stopifnot(identical(formula(glmNormal_Diff_cMMb_select), formula(glmNormal_Diff_cMMb_stepAIC)))
glmNormal_Diff_cMMb_formula <- formula(glmNormal_Diff_cMMb_select)
glmNormal_Diff_cMMb_predict <- predict(glmNormal_Diff_cMMb_select, type = "response",
                                       se = T)
glmNormal_Diff_cMMb_summary <- summary(glmNormal_Diff_cMMb_select)
glmNormal_Diff_cMMb_coeffs <- glmNormal_Diff_cMMb_summary$coefficients

save(glmNormal_Diff_cMMb_stepAIC, file = "glmNormal_Diff_cMMb_stepAIC.RData")
save(glmNormal_Diff_cMMb_formula, file = "glmNormal_Diff_cMMb_formula.RData")
save(glmNormal_Diff_cMMb_select, file = "glmNormal_Diff_cMMb_select.RData")
save(glmNormal_Diff_cMMb_predict, file = "glmNormal_Diff_cMMb_predict.RData")
save(glmNormal_Diff_cMMb_summary, file = "glmNormal_Diff_cMMb_summary.RData")
save(glmNormal_Diff_cMMb_coeffs, file = "glmNormal_Diff_cMMb_coeffs.RData")
glmNormal_Diff_cMMb_coeffs_df <- data.frame(Variable = rownames(glmNormal_Diff_cMMb_coeffs),
                                            glmNormal_Diff_cMMb_coeffs)
colnames(glmNormal_Diff_cMMb_coeffs_df) <- c("Variable", "Estimate", "Standard_error", "t-value", "Pr(>|t|)")
write.table(glmNormal_Diff_cMMb_coeffs_df,
            file = paste0(coefDir, "glmNormal_Diff_cMMb_coeffs.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

# Evaluate model goodness-of-fit based on the ratio of the
# model residual deviance to the null deviance, to give an
# R-squared value (the closer to 1 the better the model)
# https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r/46358
glmNormal_Diff_cMMb_R2 <- with(summary(glmNormal_Diff_cMMb_select), 1 - deviance/null.deviance)
print(paste0("GLM R-squared: ", glmNormal_Diff_cMMb_R2))
#GLM R-squared: 0.3941855

# Disable scientific notation for plotting
options("scipen"=100)

# Plot observed differential cM/Mb and expected differential cM/Mb by GLM (glmNormal_Diff_cMMb_select)
pdf(paste0(plotDir, "glmNormal_Diff_cMMb_observed_predicted.pdf"),
    height = 3*length(unique(dat_scaled$chr)), width = 10)
par(mar = c(5.1, 4.1, 5.1, 2.1))
par(mfrow = c(length(unique(dat_scaled$chr)), 1))
for(i in unique(dat_scaled$chr)) {
  row_indices <- which(dat_scaled$chr == i)
  dat_scaled_chr <- dat_scaled[row_indices,]
  plot(x = round((dat_scaled_chr$start+dat_scaled_chr$end)/2),
       y = glmNormal_Diff_cMMb_select$y[row_indices],
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       type = "l", lwd = 2, col = "red")
  lines(x = round((dat_scaled_chr$start+dat_scaled_chr$end)/2),
        y = glmNormal_Diff_cMMb_predict$fit[row_indices],
        lty = 2, lwd = 2, col = "black")
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = round((dat_scaled_chr$start+dat_scaled_chr$end)/2),
       labels = round(round((dat_scaled_chr$start+dat_scaled_chr$end)/2)/1e6, digits = 3))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 1, line = 2.5, cex = 1, text = paste0(i, " coordinates (Mb)"))
  mtext(side = 2, line = 2, cex = 1, text = "Differential cM/Mb")
  if(i == unique(dat_scaled$chr)[1]) {
    mtext(side = 3, line = 1.5, cex = 0.75,
          text = "Differential cM/Mb ~ IWGSC_cMMb + DMC1 + H3K4me1 + H3K27ac + H3K27me3 + H3K27me1 + CENH3 + mCG + mCHH + \nIWGSC_cMMb:DMC1 + IWGSC_cMMb:H3K4me1 + IWGSC_cMMb:H3K27me3 + DMC1:H3K4me1 + DMC1:CENH3 + H3K4me1:CENH3 + \nH3K27ac:H3K27me3 + H3K27ac:mCHH + H3K27me3:mCG + H3K27me3:mCHH + H3K27me1:mCG + CENH3:mCHH + mCG:mCHH")
    mtext(side = 3, line = 0.5, cex = 0.75,
          text = bquote(italic("R")^2 ~ "=" ~ .(round(glmNormal_Diff_cMMb_R2, digits = 2))))
  }
  box(lwd = 1.5)
  legend("topleft",
         legend = c("Observed", "Predicted"),
         text.col = c("red", "black"),
         col = "white",
         ncol = 1, cex = 1, lwd = 1.5, bty = "n")
}
dev.off()


# Partition the data into training and test sets to evaluate
# model performance on the test set

# Define training (70%) and test (30%) subsets of dat
set.seed(294502)
index <- base::sample(1:nrow(dat), size = 0.7*nrow(dat))
train <- dat[index,]
test <- dat[-index,]
dim(train)
dim(test)

# Centre and standardize the predictors
# See https://advstats.psychstat.org/book/mregression/standardization.php
# https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
x_cols <- colnames(dat)[which( !(colnames(dat) %in% c("chr", "start", "end", "Diff_cMMb") ) )]

# Apply preProcess only to the training data, and then use the results
# to scale both the training and test data
# (the predictors in the training data will be centered on 0)
pre_proc_val <- caret::preProcess(train[,x_cols], method = c("center", "scale"))

train[,x_cols] <- stats::predict(pre_proc_val, newdata = train[,x_cols])
test[,x_cols] <- stats::predict(pre_proc_val, newdata = test[,x_cols])

summary(train)

# Make model with selected predictors using the training and test sets of dat
lmNormal_Diff_cMMb_train <- lm(formula = formula(glmNormal_Diff_cMMb_stepAIC),
                                data = train)
lmNormal_Diff_cMMb_test <- lm(formula = formula(glmNormal_Diff_cMMb_stepAIC),
                              data = test)

# Define function to compute model evaluation metrics
eval_metrics <- function(model, dataFrame, predicted, observed) {
  resids <- dataFrame[,observed] - predicted
  resids2 <- resids**2
  N <- length(predicted)
  RMSE <- round(sqrt(sum(resids2, na.rm = T)/N), 2)

  r2 <- as.character(round(summary(model)$r.squared, 2))
  adj_r2 <- as.character(round(summary(model)$adj.r.squared, 2))

  print(paste0("R-squared: ", r2))
  print(paste0("Adjusted R-squared: ", adj_r2))
  print(paste0("RMSE: ", RMSE))
}

# Define function to compute R-squared and RMSE from observed and predicted values
eval_results <- function(observed, predicted, dataFrame) {
  SSE <- sum((predicted - observed)^2, na.rm = T)
  SST <- sum((observed - mean(observed))^2, na.rm = T)
  Rsquared <- 1 - SSE / SST
  RMSE <- sqrt(SSE/nrow(dataFrame))

  data.frame(Rsquared = Rsquared,
             RMSE = RMSE)
}

# Get predicted Diff_cMMb values and evaluation metrics for the model using the training set
predictions_train <- predict(lmNormal_Diff_cMMb_train, newdata = train)
eval_metrics(model = lmNormal_Diff_cMMb_train,
             dataFrame = train,
             predicted = predictions_train,
             observed = "Diff_cMMb")
#[1] "R-squared: 0.49"
#[1] "Adjusted R-squared: 0.46"
#[1] "RMSE: 115.02"
eval_results(observed = train$Diff_cMMb,
             predicted = predictions_train,
             dataFrame = train)
#   Rsquared     RMSE
#1 0.4937219 115.0204
print(paste0("Train mean Diff_cMMb: ", mean(train$Diff_cMMb, na.rm = T)))
print(paste0("Train median Diff_cMMb: ", median(train$Diff_cMMb, na.rm = T)))
print(paste0("Train range Diff_cMMb: ", range(train$Diff_cMMb, na.rm = T)))


# Get predicted Diff_cMMb values and evaluation metrics for the model using the test set
predictions_test <- predict(lmNormal_Diff_cMMb_test, newdata = test)
eval_metrics(model = lmNormal_Diff_cMMb_test,
             dataFrame = test,
             predicted = predictions_test,
             observed = "Diff_cMMb")
#[1] "R-squared: 0.66"
#[1] "Adjusted R-squared: 0.61"
#[1] "RMSE: 40.91"
eval_results(observed = test$Diff_cMMb,
             predicted = predictions_test,
             dataFrame = test)
#   Rsquared     RMSE
#1 0.6626569 40.90989
print(paste0("Test mean Diff_cMMb: ", mean(test$Diff_cMMb, na.rm = T)))
print(paste0("Test median Diff_cMMb: ", median(test$Diff_cMMb, na.rm = T)))
print(paste0("Test range Diff_cMMb: ", range(test$Diff_cMMb, na.rm = T)))


# Evaluted model performance when dat_scaled is subsetted by chromosome

pdf(paste0(plotDir, "glmNormal_Diff_cMMb_chr_observed_predicted_fit_per_chromosome.pdf"),
    height = 3*length(unique(dat_scaled$chr)), width = 10)
par(mar = c(5.1, 4.1, 5.1, 2.1))
par(mfrow = c(length(unique(dat_scaled$chr)), 1))

for(chrName in unique(dat_scaled$chr)) {

  print(chrName)
#  if( !( chrName %in% c("chr3B") ) ) {
    dat_scaled_chr <- dat_scaled[dat_scaled$chr == chrName,]
  
    glmNormal_Diff_cMMb_chr <- glm2(formula = formula(glmNormal_Diff_cMMb_stepAIC),
                                    data = dat_scaled_chr,
                                    family = gaussian(),
                                    control = glm.control(maxit = 100000))
    glmNormal_Diff_cMMb_chr_formula <- formula(glmNormal_Diff_cMMb_chr)
    glmNormal_Diff_cMMb_chr_predict <- predict(glmNormal_Diff_cMMb_chr, type = "response",
                                               se = T)
    glmNormal_Diff_cMMb_chr_summary <- summary(glmNormal_Diff_cMMb_chr)
    glmNormal_Diff_cMMb_chr_coeffs <- glmNormal_Diff_cMMb_chr_summary$coefficients
    glmNormal_Diff_cMMb_chr_R2 <- with(summary(glmNormal_Diff_cMMb_chr), 1 - deviance/null.deviance)
    print(paste0(chrName, " GLM R-squared: ", glmNormal_Diff_cMMb_chr_R2))
  
    lmNormal_Diff_cMMb_chr <- lm(formula = formula(glmNormal_Diff_cMMb_stepAIC),
                                 data = dat_scaled_chr)
    predictions_chr <- predict(lmNormal_Diff_cMMb_chr, newdata = dat_scaled_chr)
    eval_metrics(model = lmNormal_Diff_cMMb_chr,
                 dataFrame = dat_scaled_chr,
                 predicted = predictions_chr,
                 observed = "Diff_cMMb")
    print(paste0(chrName, " mean Diff_cMMb: ", mean(dat_scaled_chr$Diff_cMMb, na.rm = T)))
    print(paste0(chrName, " median Diff_cMMb: ", median(dat_scaled_chr$Diff_cMMb, na.rm = T)))
    print(paste0(chrName, " range Diff_cMMb: ",range(dat_scaled_chr$Diff_cMMb, na.rm = T)))
  
    # Plot observed differential cM/Mb and expected differential cM/Mb by GLM (glmNormal_Diff_cMMb_chr)
    plot(x = round((dat_scaled_chr$start+dat_scaled_chr$end)/2)[which(rownames(dat_scaled_chr) %in% names(glmNormal_Diff_cMMb_chr$y))],
         y = glmNormal_Diff_cMMb_chr$y,
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n",
         type = "l", lwd = 2, col = "red")
    lines(x = round((dat_scaled_chr$start+dat_scaled_chr$end)/2)[which(rownames(dat_scaled_chr) %in% names(glmNormal_Diff_cMMb_chr$y))],
          y = glmNormal_Diff_cMMb_chr_predict$fit,
          lty = 2, lwd = 2, col = "black")
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = round((dat_scaled_chr$start+dat_scaled_chr$end)/2)[which(rownames(dat_scaled_chr) %in% names(glmNormal_Diff_cMMb_chr$y))],
         labels = round(round((dat_scaled_chr$start+dat_scaled_chr$end)/2)/1e6, digits = 3)[which(rownames(dat_scaled_chr) %in% names(glmNormal_Diff_cMMb_chr$y))])
    axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
    mtext(side = 1, line = 2.5, cex = 1, text = paste0(chrName, " coordinates (Mb)"))
    mtext(side = 2, line = 2, cex = 1, text = "Differential cM/Mb")
    if(chrName == unique(dat_scaled$chr)[1]) {
      mtext(side = 3, line = 1.5, cex = 0.75,
            text = "Differential cM/Mb ~ IWGSC_cMMb + DMC1 + H3K4me1 + H3K27ac + H3K27me3 + H3K27me1 + CENH3 + mCG + mCHH + \nIWGSC_cMMb:DMC1 + IWGSC_cMMb:H3K4me1 + IWGSC_cMMb:H3K27me3 + DMC1:H3K4me1 + DMC1:CENH3 + H3K4me1:CENH3 + \nH3K27ac:H3K27me3 + H3K27ac:mCHH + H3K27me3:mCG + H3K27me3:mCHH + H3K27me1:mCG + CENH3:mCHH + mCG:mCHH")
    }
    mtext(side = 3, line = 0.5, cex = 0.75,
          text = bquote(italic("R")^2 ~ "=" ~ .(round(glmNormal_Diff_cMMb_chr_R2, digits = 2))))
    box(lwd = 1.5)
    legend("topleft",
           legend = c("Observed", "Predicted"),
           text.col = c("red", "black"),
           col = "white",
           ncol = 1, cex = 1, lwd = 1.5, bty = "n")

#  }

}

dev.off()  

#[1] "chr1A"                                                                                                                                   
#[1] "chr1A GLM R-squared: 0.811459149894994"                                                                                                  
#[1] "R-squared: 0.81"                                                                                                                         
#[1] "Adjusted R-squared: 0.78"                                                                                                                
#[1] "RMSE: 61.92"                                                                                                                             
#[1] "chr1A mean Diff_cMMb: 11.2838954317006"                                                                                                  
#[1] "chr1A median Diff_cMMb: 0"                                                                                                               
#[1] "chr1A range Diff_cMMb: -1090.19607843137"                                                                                                
#[2] "chr1A range Diff_cMMb: 1070.10710808179"                                                                                                 

#[1] "chr1B"                                                                                                                                   
#[1] "chr1B GLM R-squared: 0.830295362945433"                                                                                                  
#[1] "R-squared: 0.83"                                                                                                                         
#[1] "Adjusted R-squared: 0.78"                                                                                                                
#[1] "RMSE: 28.69"                                                                                                                             
#[1] "chr1B mean Diff_cMMb: 5.96374538272841"                                                                                                  
#[1] "chr1B median Diff_cMMb: 0"                                                                                                               
#[1] "chr1B range Diff_cMMb: -221.32763720013"                                                                                                 
#[2] "chr1B range Diff_cMMb: 609.837278106509"                                                                                                 

#[1] "chr3B"                                                                                                                                   
#[1] "chr3B GLM R-squared: 1"                                                                                                                  
#[1] "R-squared: 1"                                                                                                                            
#[1] "Adjusted R-squared: NaN"                                                                                                                 
#[1] "RMSE: 0"                                                                                                                                 
#[1] "chr3B mean Diff_cMMb: 14.3526535299572"                                                                                                  
#[1] "chr3B median Diff_cMMb: 0"                                                                                                               
#[1] "chr3B range Diff_cMMb: -34.0435953955425"                                                                                                
#[2] "chr3B range Diff_cMMb: 390.447443181818"                                                                                                 

#[1] "chr3D"                                                                                                                                   
#[1] "chr3D GLM R-squared: 1"                                                                                                                  
#[1] "R-squared: 1"                                                                                                                            
#[1] "Adjusted R-squared: NaN"
#[1] "RMSE: 0"
#[1] "chr3D mean Diff_cMMb: -52.479531348282"
#[1] "chr3D median Diff_cMMb: -0.0707662862978778"
#[1] "chr3D range Diff_cMMb: -1236.42582300128"
#[2] "chr3D range Diff_cMMb: 81.9100546892975"

#[1] "chr4A"
#[1] "chr4A GLM R-squared: 0.941898415256572"
#[1] "R-squared: 0.94"
#[1] "Adjusted R-squared: 0.84"
#[1] "RMSE: 6.99"
#[1] "chr4A mean Diff_cMMb: 6.90271389037834"
#[1] "chr4A median Diff_cMMb: 0"
#[1] "chr4A range Diff_cMMb: -23.3592020121676"
#[2] "chr4A range Diff_cMMb: 133.532140490391"

#[1] "chr5A"
#[1] "chr5A GLM R-squared: 0.992554123261895"
#[1] "R-squared: 0.99"
#[1] "Adjusted R-squared: 0.98"
#[1] "RMSE: 6.9"
#[1] "chr5A mean Diff_cMMb: -13.6323876729907"
#[1] "chr5A median Diff_cMMb: 0"
#[1] "chr5A range Diff_cMMb: -505.925250683682"
#[2] "chr5A range Diff_cMMb: 6.99293210537588"

#[1] "chr5B"
#[1] "chr5B GLM R-squared: 0.628873941603701"
#[1] "R-squared: 0.63"
#[1] "Adjusted R-squared: 0.5"
#[1] "RMSE: 33.41"
#[1] "chr5B mean Diff_cMMb: 1.60868753217813"
#[1] "chr5B median Diff_cMMb: 0"
#[1] "chr5B range Diff_cMMb: -243.859649122807"
#[2] "chr5B range Diff_cMMb: 357.236842105263"

#[1] "chr7A"
#[1] "chr7A GLM R-squared: 0.921000668332116"
#[1] "R-squared: 0.92"
#[1] "Adjusted R-squared: 0.88"
#[1] "RMSE: 57"
#[1] "chr7A mean Diff_cMMb: -12.2119877893633"
#[1] "chr7A median Diff_cMMb: 0"
#[1] "chr7A range Diff_cMMb: -941.448382126348"
#[2] "chr7A range Diff_cMMb: 883.687943262411"

#[1] "chr7B"
#[1] "chr7B GLM R-squared: 0.995307776398281"
#[1] "R-squared: 1"
#[1] "Adjusted R-squared: 0.99"
#[1] "RMSE: 16.03"
#[1] "chr7B mean Diff_cMMb: 42.7394914241486"
#[1] "chr7B median Diff_cMMb: 0.989689436654284"
#[1] "chr7B range Diff_cMMb: -287.314236485031"
#[2] "chr7B range Diff_cMMb: 1589.14342629482"
#Warning messages:
#1: In predict.lm(lmNormal_Diff_cMMb_chr, newdata = dat_scaled_chr) :
#  prediction from a rank-deficient fit may be misleading
#2: In predict.lm(lmNormal_Diff_cMMb_chr, newdata = dat_scaled_chr) :
#  prediction from a rank-deficient fit may be misleading
