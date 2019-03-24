#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./coverage_per_scaled_win_TelCen_plotAllProfiles.R both 1Mb 100 100ths 3 

#align <- "both"
#winName <- "1Mb"
#prop <- 100
#propName <- "100ths" 
#N <- 3 

args <- commandArgs(trailingOnly = T)
align <- args[1]
winName <- args[2]
prop <- as.numeric(args[3])
propName <- args[4]
N <- as.numeric(args[5])

library(parallel)

profileNames <- c(
                  #"log2_CENH3_ChIP_SRR1686799_MNase_Rep1",
                  "log2_H3K4me3_Rep1_ChIP_MNase_Rep1",
                  "log2_H3K9ac_ChIP_SRR6350667_H3_input_SRR6350669",
                  "log2_H3K4me3_ChIP_SRR6350668_H3_input_SRR6350669",
                  "log2_H3K36me3_ChIP_SRR6350670_H3_input_SRR6350669",
                  "log2_MNase_Rep1_H3_input_SRR6350669",
                  "log2_H3K27me1_Rep1_ChIP_MNase_Rep1",
                  "log2_H3K27me3_ChIP_SRR6350666_H3_input_SRR6350669",
                  "log2_H3K9me2_Rep1_ChIP_MNase_Rep1"
                 )
profileNamesPlot <- c( 
                      #"CENH3",
                      "H3K4me3",
                      "H3K9ac",
                      "H3K4me3 (IWGSC)", 
                      "H3K36me3",
                      "MNase",
                      "H3K27me1",
                      "H3K27me3",
                      "H3K9me2"
                     )
profileColours <- c(
                    #"deeppink",
                    "forestgreen",
                    "green2",
                    "grey60",
                    "darkorange2",
                    "darkcyan",
                    "firebrick1",
                    "navy",
                    "magenta3"
                   )

makeTransparent <- function(thisColour, alpha = 150)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
profileColoursTransparent <- sapply(seq_along(profileColours), function(x) {
  makeTransparent(profileColours[x])
})

# Load TelCenMatrix
TelCenDFs <- mclapply(seq_along(profileNames), function(x) {
  read.table(paste0("./", profileNames[x], "_",
             align, "_", winName, "_",
             propName, "_TelCenMatrix.txt"))
}, mc.cores = length(profileNames))

TelCenProfiles <- mclapply(seq_along(TelCenDFs), function(x) {
  as.vector(rowMeans(TelCenDFs[[x]]))
}, mc.cores = length(TelCenDFs))

TelCenSDs <- mclapply(seq_along(TelCenDFs), function(x) {
  as.vector(apply(X = TelCenDFs[[x]], MARGIN = 1, FUN = sd))
}, mc.cores = length(TelCenDFs))

# Function to plot telomere to centromere (Tel-Cen)
# log2(markChIP/markControl) profiles 
TelCenPlot <- function(xplot,
                       profiles,
                       proportions,
                       proportionsName,
                       profileColours,
                       Ylabel,
                       Ylim,
                       legendLabs,
                       legendLoc) {
  plot(xplot, profiles[[1]], type = "l", lwd = 3, col = profileColours[[1]],
       ylim = Ylim,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "Wheat")
  lines(xplot, profiles[[2]], type = "l", lwd = 3, col = profileColours[[2]])
  lines(xplot, profiles[[3]], type = "l", lwd = 3, col = profileColours[[3]])
  lines(xplot, profiles[[4]], type = "l", lwd = 3, col = profileColours[[4]])
  lines(xplot, profiles[[5]], type = "l", lwd = 3, col = profileColours[[5]])
  lines(xplot, profiles[[6]], type = "l", lwd = 3, col = profileColours[[6]])
  lines(xplot, profiles[[7]], type = "l", lwd = 3, col = profileColours[[7]])
  lines(xplot, profiles[[8]], type = "l", lwd = 3, col = profileColours[[8]])
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(1, seq(10, proportions, by = 10)),
       labels = c(expression(italic("TEL")),
                  seq(10, proportions-10, by = 10),
                  expression(italic("CEN"))))
  mtext(side = 1, line = 2.1, cex = 1,
        text = paste0("Scaled windows (", proportionsName, ")"))
  axis(side = 2,
       at = pretty(c(profiles[[1]], profiles[[2]], profiles[[3]], profiles[[4]],
                     profiles[[5]], profiles[[6]], profiles[[7]], profiles[[8]])),
       cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.1, cex = 1, text = Ylabel, col = "black")
  abline(h = 0, lwd = 1.5, lty = 1)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(profileColours),
         text.col = c(profileColours),
         text.font = c(rep(1, times = 6)),
         ncol = 2, cex = 0.7, lwd = 1.5, bty = "n")
}


# Calculate moving average of current window,
#### (N/2) previous windows (where N is even) OR
# (N/2)-0.5 previous windows (where N is odd),
# and
#### (N/2) subsequent windows (where N is even) OR
# (N/2)-0.5 subsequent windows (where N is odd)
# (the higher N is, the greater the smoothing)
stopifnot(N %% 2 != 0) 
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_TelCenProfiles <- mclapply(seq_along(TelCenProfiles), function(x) {
  filt_TelCenProfile <- stats::filter(x = TelCenProfiles[[x]],
                                      filter = f,
                                      sides = 2)
  filt_TelCenProfile[1:flank] <- filt_TelCenProfile[flank+1]
  filt_TelCenProfile[(length(filt_TelCenProfile)-flank+1):length(filt_TelCenProfile)] <- filt_TelCenProfile[(length(filt_TelCenProfile)-flank)]
  filt_TelCenProfile
}, mc.cores = length(TelCenProfiles))

pdf(paste0("./Wheat_Zscore_log2_ChIPinput_",
           align, "_", winName, "_",
           propName, "_smooth", N,
           "_TelCenProfile.pdf"),
    height = 7, width = 7)
par(mfrow = c(2, 1))
par(mar = c(3.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
TelCenPlot(xplot = 1:length(TelCenProfiles[[1]]),
           profiles = TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColoursTransparent,
           Ylabel = bquote("Log"[2]*"(ChIP/control)"),
           Ylim = c(min(c(TelCenProfiles[[1]], TelCenProfiles[[2]], TelCenProfiles[[3]], TelCenProfiles[[4]],
                          TelCenProfiles[[5]], TelCenProfiles[[6]], TelCenProfiles[[7]], TelCenProfiles[[8]])),
                    max(c(TelCenProfiles[[1]], TelCenProfiles[[2]], TelCenProfiles[[3]], TelCenProfiles[[4]],
                          TelCenProfiles[[5]], TelCenProfiles[[6]], TelCenProfiles[[7]], TelCenProfiles[[8]]))),
           legendLabs = profileNamesPlot,
           legendLoc = "top")
TelCenPlot(xplot = 1:length(filt_TelCenProfiles[[1]]),
           profiles = filt_TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColoursTransparent,
           Ylabel = bquote("Log"[2]*"(ChIP/control)"),
           Ylim = c(min(c(filt_TelCenProfiles[[1]], filt_TelCenProfiles[[2]], filt_TelCenProfiles[[3]], filt_TelCenProfiles[[4]],
                          filt_TelCenProfiles[[5]], filt_TelCenProfiles[[6]], filt_TelCenProfiles[[7]], filt_TelCenProfiles[[8]])),
                    max(c(filt_TelCenProfiles[[1]], filt_TelCenProfiles[[2]], filt_TelCenProfiles[[3]], filt_TelCenProfiles[[4]],
                          filt_TelCenProfiles[[5]], filt_TelCenProfiles[[6]], filt_TelCenProfiles[[7]], filt_TelCenProfiles[[8]]))),
           legendLabs = profileNamesPlot,
           legendLoc = "top")
dev.off()
