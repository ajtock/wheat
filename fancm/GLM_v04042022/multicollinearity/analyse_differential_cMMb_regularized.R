#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 04.04.2022

# Examine relationships between differential fancm-wild type recombination rate (cM/Mb) and
# epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq

# Build GLM with differential fancm-wild type cM/Mb as the response variable
# and epigenetic and meiotic protein signals as predictor variables

# Usage:
# ./analyse_differential_cMMb_regularized.R

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
dat <- dat[,which(colnames(dat) %in% c("diff_fancm_wt_cMMb_inter",
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
x_cols <- colnames(dat)[which(colnames(dat) != "Diff_cMMb")]
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
x_retained <- names(dat_scaled)[ which( !(names(dat_scaled) %in%
                                          c(glmNormal_Diff_cMMb_vif_gt7, "Diff_cMMb")) ) ]
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
glmNormal_Diff_cMMb_vif <- car::vif(mod = glmNormal_Diff_cMMb)
print(glmNormal_Diff_cMMb_vif)

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

# Use car::vif() to calculate generalized variance-inflation factors (GVIFs)
# for predictor variables
# https://statisticsbyjim.com/regression/multicollinearity-in-regression-analysis/
# https://www.statology.org/variance-inflation-factor-r/
# "VIFs greater than 5 represent critical levels of multicollinearity
# where the coefficients are poorly estimated, and the p-values are questionable."
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
glmNormal_Diff_cMMb_R2
#[1] 0.3941855

# Disable scientific notation for plotting
options("scipen"=100)

# Plot observed differential cM/Mb and expected differential cM/Mb by GLM (glmNormal_Diff_cMMb_select)
pdf(paste0(plotDir, "glmNormal_Diff_cMMb_observed_predicted.pdf"),
    height = 3*length(unique(dat$chr)), width = 10)
par(mfrow = c(length(unique(dat$chr)), 1))
for(i in unique(dat$chr)) {
  row_indices <- which(dat$chr == i)
  dat_chr <- dat[row_indices,]
  plot(x = round((dat_chr$start+dat_chr$end)/2),
       y = glmNormal_Diff_cMMb_select$y[row_indices],
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       type = "l", lwd = 2, col = "red")
  lines(x = round((dat_chr$start+dat_chr$end)/2),
        y = glmNormal_Diff_cMMb_predict$fit[row_indices],
        lty = 2, lwd = 2, col = "black")
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = round((dat_chr$start+dat_chr$end)/2),
       labels = round(round((dat_chr$start+dat_chr$end)/2)/1e6, digits = 3))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 1, line = 2.5, cex = 1, text = paste0(i, " coordinates (Mb)"))
  mtext(side = 2, line = 2, cex = 1, text = "Differential cM/Mb")
  if(i == unique(dat$chr)[1]) {
    mtext(side = 3, line = 2, cex = 0.75,
          text = "Differential cM/Mb ~ IWGSC_cMMb + DMC1 + H3K4me1 + H3K27ac + H3K27me3 + H3K27me1 + CENH3 + mCG + mCHH + IWGSC_cMMb:DMC1 + IWGSC_cMMb:H3K4me1 + IWGSC_cMMb:H3K27me3 + DMC1:H3K4me1 + DMC1:CENH3 + H3K4me1:CENH3 + H3K27ac:H3K27me3 + H3K27ac:mCHH + H3K27me3:mCG + H3K27me3:mCHH + H3K27me1:mCG + CENH3:mCHH + mCG:mCHH")
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
set.seed(400)
index <- base::sample(1:nrow(dat), size = 0.7*nrow(dat))
train <- dat[index,]
test <- dat[-index,]
dim(train)
dim(test)

# Centre and standardize the predictors
# See https://advstats.psychstat.org/book/mregression/standardization.php
# https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
x_cols <- colnames(dat)[which(colnames(dat) != "Diff_cMMb")]

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
eval_results(observed = train$Diff_cMMb,
             predicted = predictions_train,
             dataFrame = train)

# Get predicted Diff_cMMb values and evaluation metrics for the model using the test set
predictions_test <- predict(lmNormal_Diff_cMMb_test, newdata = test)
eval_metrics(model = lmNormal_Diff_cMMb_test,
             dataFrame = test,
             predicted = predictions_test,
             observed = "Diff_cMMb")
eval_results(observed = test$Diff_cMMb,
             predicted = predictions_test,
             dataFrame = test)







# Resources consulted on regularization to avoid model overfitting:
# https://www.geeksforgeeks.org/regularization-in-r-programming/
# https://glmnet.stanford.edu/articles/glmnet.html
# https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
# https://stackoverflow.com/questions/27580267/how-to-make-all-interactions-before-using-glmnet

# Ridge regression

# Set seed for reproducibility
set.seed(3485345)

# Get response (y) and predictor (x) variables

# Make interactions before using glmnet
# Use .*. for all interactions
f <- as.formula(y ~ .*.)

# Use model.matrix to take advantage of f
x <- model.matrix(f, dat)[,-1]

y <- dat %>% select(y) %>%
## Centre response variable (predictors will be standardized in the modelling function)
#  scale(center = TRUE, scale = F) %>%
  as.matrix()

# Perform k-fold cross-validation to select lambda
lambdas_to_try <- 10^seq(-3, 10, length.out = 100)

# Setting alpha = 0 implements ridge regression
ridge_cv <- cv.glmnet(x, y, alpha = 0,
                      lambda = lambdas_to_try,
                      standardize = TRUE,
                      nfolds = 10)

# Plot cross-validation results
# This plots the cross-validation curve (red dotted line) along with
# upper and lower standard deviation curves along the 𝜆 sequence (error bars).
# Two special values along the 𝜆 sequence are indicated by the vertical dotted lines.
# lambda.min is the value of 𝜆 that gives minimum mean cross-validated error,
# while lambda.1se is the value of 𝜆 that gives the most regularized model such that
# the cross-validated error is within one standard error of the minimum.
pdf(paste0(plotDir, "ridge_regression_10fold_cross_validation.pdf"))
plot(ridge_cv)
dev.off()

# Best cross-validated lambda
print("lambda.min")
lambda_cv <- ridge_cv$lambda.min
print(lambda_cv)

print("lamda.1se")
print(ridge_cv$lambda.1se)

# Get coefficients
coeficients_lambda.min <- coef(ridge_cv, s = "lambda.min") %>% as.matrix()
coeficients_lambda.1se <- coef(ridge_cv, s = "lambda.1se") %>% as.matrix()
identical(coeficients_lambda.min, coeficients_lambda.1se)
#[1] FALSE

# Fit final model, get its sum of squared residuals
# and multiple R-squared
model_cv <- glmnet(x, y, alpha = 0,
                   lambda = lambda_cv,
                   standardize = TRUE)
y_hat_cv <- predict(model_cv, x)
ssr_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
#ssr_cv <- sum((y - y_hat_cv) * (y - y_hat_cv))
print(ssr_cv)
rsq_ridge_cv <- cor(y, y_hat_cv)^2
print(rsq_ridge_cv)
#                  s0
#Diff_cMMb 0.09985944

# Users can explicitly control the fold that each observation is assigned to via the foldid argument.
# This is useful, for example, in using cross-validation to select a value for α
foldid <- sample(1:10, size = length(y), replace = TRUE)
cv1 <- cv.glmnet(x, y, foldid = foldid, alpha = 1)
cv.5 <- cv.glmnet(x, y, foldid = foldid, alpha = 0.5)
cv0 <- cv.glmnet(x, y, foldid = foldid, alpha = 0)

pdf(paste0(plotDir, "ridge_regression_10fold_cross_validation_alpha_comparison.pdf"))
par(mfrow = c(2,2))
plot(cv1); plot(cv.5); plot(cv0)
plot(log(cv1$lambda) , cv1$cvm , pch = 19, col = "red",
xlab = "log(Lambda)", ylab = cv1$name)
points(log(cv.5$lambda), cv.5$cvm, pch = 19, col = "grey")
points(log(cv0$lambda) , cv0$cvm , pch = 19, col = "blue")
legend("topleft", legend = c("alpha= 1", "alpha= .5", "alpha= 0"),
pch = 19, col = c("red","grey","blue"))
dev.off()


# Lasso (least absolute shrinkage and selection operator) regression

# Set seed for reproducibility
set.seed(3485345)

# Get response (y) and predictor (x) variables

# Make interactions before using glmnet
# Use .*. for all interactions
f <- as.formula(y ~ .*.)

# Use model.matrix to take advantage of f
x <- model.matrix(f, dat)[,-1]

y <- dat %>% select(y) %>%
## Centre response variable (predictors will be standardized in the modelling function)
#  scale(center = TRUE, scale = F) %>%
  as.matrix()

# Perform k-fold cross-validation to select lambda
lambdas_to_try <- 10^seq(-3, 10, length.out = 100)

# Setting alpha = 1 implements lasso regression
lasso_cv <- cv.glmnet(x, y, alpha = 1,
                      lambda = lambdas_to_try,
                      standardize = TRUE,
                      nfolds = 10)

# Plot cross-validation results
# This plots the cross-validation curve (red dotted line) along with
# upper and lower standard deviation curves along the 𝜆 sequence (error bars).
# Two special values along the 𝜆 sequence are indicated by the vertical dotted lines.
# lambda.min is the value of 𝜆 that gives minimum mean cross-validated error,
# while lambda.1se is the value of 𝜆 that gives the most regularized model such that
# the cross-validated error is within one standard error of the minimum.
pdf(paste0(plotDir, "lasso_regression_10fold_cross_validation.pdf"))
plot(lasso_cv)
dev.off()

# Best cross-validated lambda
print("lambda.min")
lambda_cv <- lasso_cv$lambda.min
print(lambda_cv)

print("lamda.1se")
print(lasso_cv$lambda.1se)

# Get coefficients
coeficients_lambda.min <- coef(lasso_cv, s = "lambda.min") %>% as.matrix()
coeficients_lambda.1se <- coef(lasso_cv, s = "lambda.1se") %>% as.matrix()
identical(coeficients_lambda.min, coeficients_lambda.1se)
#[1] TRUE

# Fit final model, get its sum of squared residuals
# and multiple R-squared
model_cv <- glmnet(x, y, alpha = 1,
                   lambda = lambda_cv,
                   standardize = TRUE)
y_hat_cv <- predict(model_cv, x)
ssr_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
#ssr_cv <- sum((y - y_hat_cv) * (y - y_hat_cv))
print(ssr_cv)
rsq_lasso_cv <- cor(y, y_hat_cv)^2
# Warning message:
# In cor(y, y_hat_cv) : the standard deviation is zero
print(rsq_lasso_cv)
#  s0
#y NA

# Elastic net regression
# Elastic net regression can be stated as the convex combination of the lasso and ridge regression.
# Can be done with glmnet::glmnet() or caret::train()

## glmnet

# Set seed for reproducibility
set.seed(3485345)

# Get response (y) and predictor (x) variables

# Make interactions before using glmnet
# Use .*. for all interactions
f <- as.formula(y ~ .*.)

# Use model.matrix to take advantage of f
x <- model.matrix(f, dat)[,-1]

y <- dat %>% select(y) %>%
## Centre response variable (predictors will be standardized in the modelling function)
#  scale(center = TRUE, scale = F) %>%
  as.matrix()

# Perform k-fold cross-validation to select lambda
lambdas_to_try <- 10^seq(-3, 10, length.out = 100)

# Setting alpha = 0 implements ridge regression
enet_cv <- cv.glmnet(x, y, alpha = 0.5,
                     lambda = lambdas_to_try,
                     standardize = TRUE,
                     nfolds = 10)

# Plot cross-validation results
# This plots the cross-validation curve (red dotted line) along with
# upper and lower standard deviation curves along the 𝜆 sequence (error bars).
# Two special values along the 𝜆 sequence are indicated by the vertical dotted lines.
# lambda.min is the value of 𝜆 that gives minimum mean cross-validated error,
# while lambda.1se is the value of 𝜆 that gives the most regularized model such that
# the cross-validated error is within one standard error of the minimum.
pdf(paste0(plotDir, "enet_regression_10fold_cross_validation.pdf"))
plot(enet_cv)
dev.off()

# Best cross-validated lambda
print("lambda.min")
lambda_cv <- enet_cv$lambda.min
print(lambda_cv)

print("lamda.1se")
print(enet_cv$lambda.1se)

# Get coefficients
coeficients_lambda.min <- coef(enet_cv, s = "lambda.min") %>% as.matrix()
coeficients_lambda.1se <- coef(enet_cv, s = "lambda.1se") %>% as.matrix()
identical(coeficients_lambda.min, coeficients_lambda.1se)
#[1] FALSE

# Fit final model, get its sum of squared residuals
# and multiple R-squared
model_cv <- glmnet(x, y,
                   lambda = lambda_cv,
                   standardize = TRUE)
y_hat_cv <- predict(model_cv, x)
ssr_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
#ssr_cv <- sum((y - y_hat_cv) * (y - y_hat_cv))
print(ssr_cv)
rsq_ridge_cv <- cor(y, y_hat_cv)^2
print(rsq_ridge_cv)
#                  s0
#Diff_cMMb 0.09985944




## caret
# Set training control
train_control <- trainControl(method = "repeatedcv",
                              number = 10,
                              repeats = 5,
                              search = "random",
                              verboseIter = TRUE)
  
# Train the model
elastic_net_model <- train(y ~ .,
                           data = cbind(y, x),
                           method = "glmnet",
                           preProcess = c("center", "scale"),
                           tuneLength = 10,
                           trControl = train_control)


# Evaluation metrics function
eval_metrics = function(model, df, predictions, target) {
  resids = df[,target] - predictions
  resids2 = resids**2
  N = length(predictions)
  r2 = as.character(round(summary(model)$r.squared, 2))
  adj_r2 = as.character(round(summary(model)$adj.r.squared, 2))
  print(adj_r2) #Adjusted R-squared
  print(as.character(round(sqrt(sum(resids2)/N), 2))) #RMSE
}

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  Rsquared <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))

  # Model performance metrics
  data.frame(RMSE = RMSE,
             Rsquared = Rsquared)
}

eval_results(true = y,
             predicted = predict(elastic_net_model, x),
             df = dat)

# Check multiple R-squared
y_hat_enet <- predict(elastic_net_model, x)
rsq_enet <- cor(y, y_hat_enet)^2
 
print(y_hat_enet)
print(rsq_enet)

predict(elastic_net_model$finalModel, type = "coefficients", s = elastic_net_model$bestTune$lambda) 

