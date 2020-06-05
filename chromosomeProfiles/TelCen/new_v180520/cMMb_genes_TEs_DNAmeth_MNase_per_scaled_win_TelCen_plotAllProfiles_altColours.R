#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./cMMb_genes_TEs_DNAmeth_MNase_per_scaled_win_TelCen_plotAllProfiles_altColours.R 1Mb 100 100ths 3 

#winName <- "1Mb"
#prop <- 100
#propName <- "100ths" 
#N <- 3 

args <- commandArgs(trailingOnly = T)
winName <- args[1]
prop <- as.numeric(args[2])
propName <- args[3]
N <- as.numeric(args[4])

library(parallel)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

profileNames <- c(
                  "cMMb_1Mb",
                  "genes_1Mb",
                  "TEs_1Mb_superfamily_Mariner_DTT",
                  "TEs_1Mb_superfamily_Gypsy_LTR_RLG",
                  "BSseq_Rep8a_SRR6792678_mCG_1Mb",
                  "BSseq_Rep8a_SRR6792678_mCHG_1Mb",
                  "BSseq_Rep8a_SRR6792678_mCHH_1Mb",
                  "MNase_Rep1_both_1Mb"
                 )
profileNamesPlot <- c( 
                      "cM/Mb",
                      "Genes",
                      "Mariner TEs",
                      "Gypsy LTR TEs",
                      "mCG",
                      "mCHG",
                      "mCHH",
                      "MNase"
                     )
profileColours <- c(
                    "red",
                    "red4",
                    "magenta3",
                    "purple4",
                    "navy",
                    "blue",
                    "deepskyblue1",
                    "darkcyan"
                   )
profileNames <- rev(profileNames)
profileNamesPlot <- rev(profileNamesPlot)
profileColours <- rev(profileColours)

makeTransparent <- function(thisColour, alpha = 210)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
profileColours <- sapply(seq_along(profileColours), function(x) {
  makeTransparent(profileColours[x])
})

# Load TelCenMatrix
TelCenDFs <- mclapply(seq_along(profileNames), function(x) {
  read.table(paste0(profileNames[x], "_",
                    propName, "_TelCenMatrix.txt"))
}, mc.cores = length(profileNames))

TelCenProfiles <- mclapply(seq_along(TelCenDFs), function(x) {
  as.vector(rowMeans(TelCenDFs[[x]], na.rm = T))
}, mc.cores = length(TelCenDFs))

TelCenSDs <- mclapply(seq_along(TelCenDFs), function(x) {
  as.vector(apply(X = TelCenDFs[[x]], MARGIN = 1, FUN = sd, na.rm = T))
}, mc.cores = length(TelCenDFs))

# Function to plot telomere to centromere (Tel-Cen)
# log2(markChIP/markControl) profiles 
TelCenPlot <- function(xplot,
                       profiles,
                       proportions,
                       proportionsName,
                       profileColours,
                       profileColoursLeg,
                       Ylabel,
                       legendLabs,
                       legendLocx,
                       legendLocy) {
  # Right y-axis
  plot(xplot, profiles[[1]], type = "l", lwd = 3, lty = 5, col = profileColours[1],
       ylim = c(min(profiles[[1]], na.rm = T),
                max(profiles[[1]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[1], col.axis = profileColours[1], line = 0.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[2]], type = "l", lwd = 3, lty = 1, col = profileColours[2],
       ylim = c(min(profiles[[2]], na.rm = T),
                max(profiles[[2]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[2], col.axis = profileColours[2], line = 2.0)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[3]], type = "l", lwd = 3, lty = 5, col = profileColours[3],
       ylim = c(min(profiles[[3]], na.rm = T),
                max(profiles[[3]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[3], col.axis = profileColours[3], line = 3.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[4]], type = "l", lwd = 3, lty = 1, col = profileColours[4],
       ylim = c(min(profiles[[4]], na.rm = T),
                max(profiles[[4]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[4], col.axis = profileColours[4], line = 5.0)

  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.0, adj = c(0.5, -10.5), labels = Ylabel, xpd = NA, srt = -90, col = "black")

  # Left y-axis
  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[5]], type = "l", lwd = 3, lty = 5, col = profileColours[5],
       ylim = c(min(profiles[[5]], na.rm = T),
                max(profiles[[5]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[5], col.axis = profileColours[5], line = 0.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[6]], type = "l", lwd = 3, lty = 1, col = profileColours[6],
       ylim = c(min(profiles[[6]], na.rm = T),
                max(profiles[[6]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[6], col.axis = profileColours[6], line = 2.0)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[7]], type = "l", lwd = 3, lty = 5, col = profileColours[7],
       ylim = c(min(profiles[[7]], na.rm = T),
                max(profiles[[7]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[7], col.axis = profileColours[7], line = 3.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[8]], type = "l", lwd = 3, lty = 1, col = profileColours[8],
       ylim = c(min(profiles[[8]], na.rm = T),
                max(profiles[[8]], na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0, col = profileColours[8], col.axis = profileColours[8], line = 5.0)

  mtext(side = 2, line = 8, cex = 1, text = Ylabel, col = "black")

  par(mgp = c(3, 0.5, 0))
  axis(side = 1, cex.axis = 1, lwd = 2.0, lwd.tick = 2.0,
       at = c(1, seq(10, proportions, by = 10)),
       labels = c(expression(italic("TEL")),
                  seq(10, proportions-10, by = 10),
                  expression(italic("CEN"))))
  mtext(side = 1, line = 1.5, cex = 1,
        text = paste0("Scaled windows (", proportionsName, ")"))

  box(lwd = 2.0)
  legend(legendLocx, legendLocy,
         legend = legendLabs,
         col = c("white"),
         text.col = c(profileColoursLeg),
         text.font = c(rep(1, times = 9)),
         ncol = 2, cex = 0.6, lwd = 1.5, bty = "n")
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

pdf(paste0(plotDir, "cMMb_genes_TEs_DNAmeth_MNase_",
           winName, "_",
           propName, "_smooth", N,
           "_TelCenProfile_v220520_altColours.pdf"),
    height = 7, width = 12)
par(mfrow = c(2, 1))
par(mar = c(3.1, 10.1, 2.1, 10.1))
par(mgp = c(3, 0.5, 0))
TelCenPlot(xplot = 1:length(TelCenProfiles[[1]]),
           profiles = TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColours,
           profileColoursLeg = rev(profileColours),
           Ylabel = bquote(""),
           legendLabs = rev(profileNamesPlot),
           legendLocx = 8, legendLocy = 2.7)
TelCenPlot(xplot = 1:length(filt_TelCenProfiles[[1]]),
           profiles = filt_TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColours,
           profileColoursLeg = rev(profileColours),
           Ylabel = bquote(""),
           legendLabs = rev(profileNamesPlot),
           legendLocx = 8, legendLocy = 2.1)
dev.off()
