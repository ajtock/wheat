#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./coverage_per_scaled_win_TelCen_arabidopsis_plot2Profiles.R log2_WT_H3K27me1_ChIP_WT_histone_input H3K27me1 log2_WT_H3K27me3_ChIP_WT_histone_input H3K27me3 both 10kb 100 100ths 3 firebrick1 navy

args <- commandArgs(trailingOnly = T)
profileNameA <- args[1]
profileNamePlotA <- args[2]
profileNameB <- args[3]
profileNamePlotB <- args[4]
align <- args[5]
winName <- args[6]
prop <- as.numeric(args[7])
propName <- args[8]
N <- as.numeric(args[9])
colourA <- args[10]
colourB <- args[11]

library(parallel)

profileNames <- c(
                  profileNameA,
                  profileNameB
                 )
profileNamesPlot <- c( 
                      profileNamePlotA,
                      profileNamePlotB
                     )
profileColours <- c(
                    colourA,
                    colourB
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
       main = bquote(atop(bold("Arabidopsis"),
                          italic("r"[s]) ~ " = " ~ .(round(cor(profiles[[1]],
                                                               profiles[[2]],
                                                               method = "spearman"),
                                                           digits = 2)))))
  lines(xplot, profiles[[2]], type = "l", lwd = 3, col = profileColours[[2]])
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(1, seq(10, proportions, by = 10)),
       labels = c(expression(italic("TEL")),
                  seq(10, proportions-10, by = 10),
                  expression(italic("CEN"))))
  mtext(side = 1, line = 2.1, cex = 1,
        text = paste0("Scaled windows (", proportionsName, ")"))
  axis(side = 2,
       at = pretty(c(profiles[[1]], profiles[[2]])),
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

pdf(paste0("./Arabidopsis_Zscore_",
           profileNameA, "_", profileNameB, "_",
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
           Ylim = c(min(c(TelCenProfiles[[1]], TelCenProfiles[[2]])),
                    max(c(TelCenProfiles[[1]], TelCenProfiles[[2]]))),
           legendLabs = profileNamesPlot,
           legendLoc = "top")
TelCenPlot(xplot = 1:length(filt_TelCenProfiles[[1]]),
           profiles = filt_TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColoursTransparent,
           Ylabel = bquote("Log"[2]*"(ChIP/control)"),
           Ylim = c(min(c(filt_TelCenProfiles[[1]], filt_TelCenProfiles[[2]])),
                    max(c(filt_TelCenProfiles[[1]], filt_TelCenProfiles[[2]]))),
           legendLabs = profileNamesPlot,
           legendLoc = "top")
dev.off()
