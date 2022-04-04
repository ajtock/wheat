#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 04.04.2022

# Examine relationships between differential fancm-wild type recombination rate (cM/Mb) and
# epigenetic marks and meiotic proteins, detected by ChIP-seq and BS-seq

# Build GLM with differential fancm-wild type cM/Mb as the response variable
# and epigenetic and meiotic protein signals as predictor variables

# Usage:
# ./analyse_differential_incl_interactions.R

options(stringsAsFactors = F)
library(glmnet)
library(caret)
library(dplyr)
#library(fitdistrplus) # descdist, plotdist, fitdist included
#library(glm2)
#library(MASS) # glm.nb included; MASS is also loaded as a dependency of fitdistrplus
#library(pscl) # zeroinfl included
#library(vcd) # goodfit included
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

# Load cM/Mb data
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
                                       "CENH3", "mCG", "mCHG", "mCHH"))]
#                                       "CENH3", "mCG", "mCHG", "mCHH", "cMMb_10Mb"))]

colnames(dat)[which(colnames(dat) == "diff_fancm_wt_cMMb_inter")] <- "y"
#colnames(dat)[which(colnames(dat) == "cMMb_10Mb")] <- "IWGSC_cMMb"


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

