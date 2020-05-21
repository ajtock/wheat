#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./log2ChIPcontrol_coverage_per_scaled_win_TelCen_plotAllProfiles_noZscore.R both 1Mb 100 100ths 3 

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

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

profileNames <- c(
                  "log2_DMC1_Rep1_ChIP_H3_input_SRR6350669",
                  "log2_ASY1_CS_Rep1_ChIP_H3_input_SRR6350669",
                  "log2_H3K4me1_Rep1_ChIP_SRR8126618_H3_input_SRR6350669",
                  "log2_H3K4me3_Rep1_ChIP_MNase_Rep1",
                  "log2_H3K27ac_Rep1_ChIP_SRR8126621_H3_input_SRR6350669",
                  "log2_H3K27me3_ChIP_SRR6350666_H3_input_SRR6350669",
                  "log2_H3K36me3_ChIP_SRR6350670_H3_input_SRR6350669",
                  "log2_H3K27me1_Rep1_ChIP_MNase_Rep1",
                  "log2_H3K9me2_Rep1_ChIP_MNase_Rep1",
                  "log2_CENH3_ChIP_SRR1686799_H3_input_SRR6350669"
                 )
profileNamesPlot <- c( 
                      "DMC1",
                      "ASY1",
                      "H3K4me1",
                      "H3K4me3",
                      "H3K27ac",
                      "H3K27me3",
                      "H3K36me3",
                      "H3K27me1",
                      "H3K9me2",
                      "CENH3"
                     )
profileColours <- c(
                    "green2",
                    "purple4",
                    "goldenrod1",
                    "forestgreen",
                    "darkorange2",
                    "dodgerblue",
                    "blue",
                    "firebrick1",
                    "magenta3",
                    "navy"
                   )
profileNames <- rev(profileNames)
profileNamesPlot <- rev(profileNamesPlot)
profileColours <- rev(profileColours)

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
  read.table(paste0(profileNames[x], "_",
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
                       profileColoursLeg,
                       Ylabel,
                       legendLabs,
                       legendLoc) {
  # Right y-axis
  plot(xplot, profiles[[1]], type = "l", lwd = 3, col = profileColours[1],
       ylim = c(-max(c(min(profiles[[1]], na.rm = T)*-1, max(profiles[[1]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[1]], na.rm = T)*-1, max(profiles[[1]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5, col = profileColours[1], col.axis = profileColours[1], line = 0.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[2]], type = "l", lwd = 3, col = profileColours[2],
       ylim = c(-max(c(min(profiles[[2]], na.rm = T)*-1, max(profiles[[2]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[2]], na.rm = T)*-1, max(profiles[[2]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5, col = profileColours[2], col.axis = profileColours[2], line = 2.0)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[3]], type = "l", lwd = 3, col = profileColours[3],
       ylim = c(-max(c(min(profiles[[3]], na.rm = T)*-1, max(profiles[[3]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[3]], na.rm = T)*-1, max(profiles[[3]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5, col = profileColours[3], col.axis = profileColours[3], line = 3.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[4]], type = "l", lwd = 3, col = profileColours[4],
       ylim = c(-max(c(min(profiles[[4]], na.rm = T)*-1, max(profiles[[4]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[4]], na.rm = T)*-1, max(profiles[[4]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5, col = profileColours[4], col.axis = profileColours[4], line = 5.0)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[5]], type = "l", lwd = 3, col = profileColours[5],
       ylim = c(-max(c(min(profiles[[5]], na.rm = T)*-1, max(profiles[[5]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[5]], na.rm = T)*-1, max(profiles[[5]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5, col = profileColours[5], col.axis = profileColours[5], line = 6.5)

  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.0, adj = c(0.5, -10.5), labels = Ylabel, xpd = NA, srt = -90, col = "black")

  # Left y-axis
  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[6]], type = "l", lwd = 3, col = profileColours[6],
       ylim = c(-max(c(min(profiles[[6]], na.rm = T)*-1, max(profiles[[6]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[6]], na.rm = T)*-1, max(profiles[[6]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5, col = profileColours[6], col.axis = profileColours[6], line = 0.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[7]], type = "l", lwd = 3, col = profileColours[7],
       ylim = c(-max(c(min(profiles[[7]], na.rm = T)*-1, max(profiles[[7]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[7]], na.rm = T)*-1, max(profiles[[7]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5, col = profileColours[7], col.axis = profileColours[7], line = 2.0)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[8]], type = "l", lwd = 3, col = profileColours[8],
       ylim = c(-max(c(min(profiles[[8]], na.rm = T)*-1, max(profiles[[8]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[8]], na.rm = T)*-1, max(profiles[[8]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5, col = profileColours[8], col.axis = profileColours[8], line = 3.5)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[9]], type = "l", lwd = 3, col = profileColours[9],
       ylim = c(-max(c(min(profiles[[9]], na.rm = T)*-1, max(profiles[[9]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[9]], na.rm = T)*-1, max(profiles[[9]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5, col = profileColours[9], col.axis = profileColours[9], line = 5.0)

  par(new = T, mgp = c(3, 0.5, 0))
  plot(xplot, profiles[[10]], type = "l", lwd = 3, col = profileColours[10],
       ylim = c(-max(c(min(profiles[[10]], na.rm = T)*-1, max(profiles[[10]], na.rm = T)), na.rm = T),
                max(c(min(profiles[[10]], na.rm = T)*-1, max(profiles[[10]], na.rm = T)), na.rm = T)),
       xlab = "", ylab = "", main = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5, col = profileColours[10], col.axis = profileColours[10], line = 6.5)

  mtext(side = 2, line = 8, cex = 1, text = Ylabel, col = "black")

  par(mgp = c(3, 0.5, 0))
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(1, seq(10, proportions, by = 10)),
       labels = c(expression(italic("TEL")),
                  seq(10, proportions-10, by = 10),
                  expression(italic("CEN"))))
  mtext(side = 1, line = 1.5, cex = 1,
        text = paste0("Scaled windows (", proportionsName, ")"))

  abline(h = 0, lwd = 1.5, lty = 1)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c("white"),
         text.col = c(profileColoursLeg),
         text.font = c(rep(1, times = 9)),
         ncol = 2, cex = 0.8, lwd = 1.5, bty = "n")
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

pdf(paste0(plotDir, "log2ChIPcontrol_",
           align, "_", winName, "_",
           propName, "_smooth", N,
           "_TelCenProfile_v210520_reordered.pdf"),
    height = 7, width = 12)
par(mfrow = c(2, 1))
par(mar = c(3.1, 10.1, 2.1, 10.1))
par(mgp = c(3, 0.5, 0))
TelCenPlot(xplot = 1:length(TelCenProfiles[[1]]),
           profiles = TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColoursTransparent,
           profileColoursLeg = rev(profileColoursTransparent),
           Ylabel = bquote("Log"[2]*"(ChIP/control)"),
           legendLabs = rev(profileNamesPlot),
           legendLoc = "top")
TelCenPlot(xplot = 1:length(filt_TelCenProfiles[[1]]),
           profiles = filt_TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColoursTransparent,
           profileColoursLeg = rev(profileColoursTransparent),
           Ylabel = bquote("Log"[2]*"(ChIP/control)"),
           legendLabs = rev(profileNamesPlot),
           legendLoc = "top")
dev.off()
