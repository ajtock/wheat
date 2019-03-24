#!/applications/R/R-3.5.0/bin/Rscript

# Plot library-size-normalized average coverage profiles around features and
# around equivalent random loci

# Usage:
# ./featureProfilesPlot_ChIPseq_MNaseSeq_log2ChIPinput.R both hotspot_chr3B_9188023_9189966 20bp 20 10kb 10000 280219

align <- "both"
featureName <- "hotspot_chr3B_9188023_9189966"
binName <- "20bp"
binSize <- 20
flankName <- "10kb"
flankSize <- 10000
date <- "280219"

args <- commandArgs(trailingOnly = T)
align <- args[1]
featureName <- args[2]
binName <- args[3]
binSize <- as.numeric(args[4])
flankName <- args[5]
flankSize <- as.numeric(args[6])
date <- as.character(args[7])

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(parallel)

inDir <- "./matrices/"
plotDir <- "./plots/"

# Sample names and directories
profileNames <- c(
                  "H3K4me3_Rep1_ChIP",
                  "H3K9ac_ChIP_SRR6350667",
                  "H3K4me3_ChIP_SRR6350668",
                  "H3K36me3_ChIP_SRR6350670",
                  "MNase_Rep1",
                  "H3K27me1_Rep1_ChIP",
                  "H3K27me3_ChIP_SRR6350666",
                  "H3K9me2_Rep1_ChIP",
                  "H3_input_SRR6350669"
                 )

profileNamesPlot <- c(
                      "H3K4me3",
                      "H3K9ac",
                      "H3K4me3 (IWGSC)",
                      "H3K36me3",
                      "MNase",
                      "H3K27me1",
                      "H3K27me3",
                      "H3K9me2",
                      "Input"
                     )
profileColours <- c(
                    "forestgreen",
                    "green2",
                    "grey60",
                    "darkorange2",
                    "darkcyan",
                    "firebrick1",
                    "navy",
                    "magenta3",
                    "black"
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

# Load coverage matrices and mean of each column
featureColMeanslist <- mclapply(seq_along(profileNames), function(x) {
  as.vector( colMeans( 
                      read.table(paste0(inDir,
                                        profileNames[x],
                                        "_MappedOn_wheat_v1.0_lowXM_",
                                        align, "_sort_norm_",
                                        featureName, "_matrix_bin", binName,
                                        "_flank", flankName, ".tab"),
                                 skip = 3)
                      , na.rm = T) )
}, mc.cores = length(profileNames))

profiles <- list(
                 log2((featureColMeanslist[[1]]+1)/(featureColMeanslist[[5]]+1)),
                 log2((featureColMeanslist[[2]]+1)/(featureColMeanslist[[9]]+1)),
                 log2((featureColMeanslist[[3]]+1)/(featureColMeanslist[[9]]+1)),
                 log2((featureColMeanslist[[4]]+1)/(featureColMeanslist[[9]]+1)),
                 log2((featureColMeanslist[[5]]+1)/(featureColMeanslist[[9]]+1)),
                 log2((featureColMeanslist[[6]]+1)/(featureColMeanslist[[5]]+1)),
                 log2((featureColMeanslist[[7]]+1)/(featureColMeanslist[[9]]+1)),
                 log2((featureColMeanslist[[8]]+1)/(featureColMeanslist[[5]]+1))
                )

# Function to plot coverage profile of multiple chromatin marks around feature
plotFeatureCov <- function(xplot,
                           profiles,
                           profileColours,
                           Ylabel,
                           Ylim,
                           flankSize, binSize,
                           flankLabL, flankLabR,
                           featureStartLab, featureEndLab,
                           legendLabs,
                           legendLoc) {
  # Feature loci
  plot(xplot, profiles[[1]], type = "l", lwd = 3, col = profileColours[[1]],
       ylim = Ylim,
       xaxt = "n", ann = F)
  lines(xplot, profiles[[2]], type = "l", lwd = 3, col = profileColours[[2]])
  lines(xplot, profiles[[3]], type = "l", lwd = 3, col = profileColours[[3]])
  lines(xplot, profiles[[4]], type = "l", lwd = 3, col = profileColours[[4]])
  lines(xplot, profiles[[5]], type = "l", lwd = 3, col = profileColours[[5]])
  lines(xplot, profiles[[6]], type = "l", lwd = 3, col = profileColours[[6]])
  lines(xplot, profiles[[7]], type = "l", lwd = 3, col = profileColours[[7]])
  lines(xplot, profiles[[8]], type = "l", lwd = 3, col = profileColours[[8]])
  axis(side = 2,
       at = pretty(c(profiles[[1]], profiles[[2]], profiles[[3]], profiles[[4]],
                     profiles[[5]], profiles[[6]], profiles[[7]], profiles[[8]])),
       cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 1, text = Ylabel, col = "black")
  axis(side = 1, labels = c("", "", "", ""), cex.axis = 1, lwd.tick = 1.5,
       at = c(1,
              (flankSize/binSize)+1,
              length(profiles[[1]])-(flankSize/binSize),
              length(profiles[[1]])))
  mtext(side = 1, line = 1, cex = 0.7,
        text = c(flankLabL,
                 featureStartLab,
                 featureEndLab,
                 flankLabR),
        at = c(1,
               (flankSize/binSize)+1,
               length(profiles[[1]])-(flankSize/binSize),
               length(profiles[[1]])))
  abline(v = c((flankSize/binSize)+1,
               length(profiles[[1]])-(flankSize/binSize)), lty = 3, lwd = 2)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(profileColours),
         text.col = c(profileColours),
         text.font = c(rep(1, times = 8)),
         ncol = 2, cex = 0.7, lwd = 1.5, bty = "n")
}

# Plot
pdf(paste0(plotDir, "Wheat_", featureName, "_profiles_",
           align, "_bin", binName, "_flank", flankName,
           "_log2ChIPinput_v", date, ".pdf"),
   height = 2.5, width = 10)
par(mfrow = c(1, 1))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

plotFeatureCov(xplot = seq_along(profiles[[1]]),
               profiles = profiles,
               profileColours = profileColoursTransparent,
               Ylabel = bquote("Log"[2]*"(ChIP/control)"),
               Ylim = c(min(c(profiles[[1]], profiles[[2]], profiles[[3]], profiles[[4]],
                              profiles[[5]], profiles[[6]], profiles[[7]], profiles[[8]])),
                        max(c(profiles[[1]], profiles[[2]], profiles[[3]], profiles[[4]],
                              profiles[[5]], profiles[[6]], profiles[[7]], profiles[[8]]))),
               flankSize = flankSize, binSize = binSize,
               flankLabL = "-10 kb", flankLabR = "+10 kb",
               featureStartLab = featureStartLab, featureEndLab = featureEndLab,
               legendLabs = profileNamesPlot[1:8],
               legendLoc = "top")
dev.off()
