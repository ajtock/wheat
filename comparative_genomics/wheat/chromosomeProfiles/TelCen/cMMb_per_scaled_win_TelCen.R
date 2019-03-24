#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./cMMb_per_scaled_win_TelCen.R 200 1Mb 1000000 100 100ths red

minMarkerDist <- 200
winName <- "1Mb"
winSize <- 1000000
prop <- 100
propName <- "100ths" 
profileColour <- "red"

args <- commandArgs(trailingOnly = T)
minMarkerDist <- as.numeric(args[1])
winName <- args[2]
winSize <- as.numeric(args[3])
prop <- as.numeric(args[4])
propName <- args[5]
profileColour <- args[6]

library(parallel)
library(segmentSeq)

makeTransparent <- function(thisColour, alpha = 150)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
profileColour <- sapply(seq_along(profileColour), function(x) {
  makeTransparent(profileColour[x])
})

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
centromere <- centromereStart+((centromereEnd-centromereStart)/2)


## cMMb profile
inDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/"

cMMb <- read.table(paste0(inDir,
                          "cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                          minMarkerDist, "bp_", winName, ".txt"),
                   header = T)
chrProfiles <- mclapply(seq_along(chrs), function(x) {
  cMMb[cMMb$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate average cMMb in proportionally scaled
# windows along left and right chromosome arms
TelCenMatrix <- NULL
for(x in seq_along(chrs)) {
  print(chrs[x])
  # Left arm
  seqWindowStartsL <- as.integer(seq(from = 1,
                                     to = centromere[x],
                                     by = centromere[x]/prop))
  # Check that last window start coordinate is
  # < centromere[x]-1
  stopifnot(seqWindowStartsL[length(seqWindowStartsL)] < centromere[x]-1)
  seqWindowEndsL <- c(seqWindowStartsL[2:length(seqWindowStartsL)]-1,
                      centromere[x]-1)
  windowGRangesL <- GRanges(seqnames = chrs[x],
                            ranges = IRanges(start = seqWindowStartsL,
                                             end = seqWindowEndsL),
                            strand = "*")
  print(windowGRangesL)
  # Check that last window end coordinate == centromere[x]-1
  # (i.e., most proximal coordinate)
  stopifnot(end(windowGRangesL[length(windowGRangesL)]) == centromere[x]-1)

  # Right arm
  seqWindowStartsR <- as.integer(seq(from = centromere[x],
                                     to = chrLens[x],
                                     by = ((chrLens[x]-centromere[x])+1)/prop))
  # Check that last window start coordinate is
  # < chrLens[x]
  stopifnot(seqWindowStartsR[length(seqWindowStartsR)] < chrLens[x])
  seqWindowEndsR <- c(seqWindowStartsR[2:length(seqWindowStartsR)]-1,
                      chrLens[x])
  # Reverse window order so distal windows come first
  windowGRangesR <- rev(GRanges(seqnames = chrs[x],
                                 ranges = IRanges(start = seqWindowStartsR,
                                                  end = seqWindowEndsR),
                                 strand = "*"))
  print(windowGRangesR)
  # Check that first window end coordinate == chrLens[x]
  # (i.e., most distal coordinate)
  stopifnot(end(windowGRangesR[1]) == chrLens[x])

  # Create GRanges object of winName windows
  # corresponding to winName windowed cMMb values
  cMMbWindowStarts <- chrProfiles[[x]]$window
  cMMbWindowEnds <- c((chrProfiles[[x]]$window[2:length(chrProfiles[[x]]$window)])-1,
                       chrLens[x])
  cMMbWindowGRanges <- GRanges(seqnames = chrs[x],
                               ranges = IRanges(start = cMMbWindowStarts,
                                                end = cMMbWindowEnds),
                               strand = "*")

  # Calculate mean cMMb in each
  # scaled window using cMMb values in winName windows
  overlapsL <- getOverlaps(coordinates = windowGRangesL,
                           segments = cMMbWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  overlapsR <- getOverlaps(coordinates = windowGRangesR,
                           segments = cMMbWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  scaledWinAvgCovL <- sapply(overlapsL, function(y) {
                        mean(chrProfiles[[x]]$cMMb[y], na.rm = T)
                      })
  scaledWinAvgCovR <- sapply(overlapsR, function(y) {
                        mean(chrProfiles[[x]]$cMMb[y], na.rm = T)
                      })
  scaledWinAvgCovLR <- cbind(scaledWinAvgCovL,
                             scaledWinAvgCovR)
  TelCenMatrix <- cbind(TelCenMatrix, scaledWinAvgCovLR)
}
write.table(TelCenMatrix,
            file = paste0("./cMMb_iwgsc_refseqv1.0_minInterMarkerDist",
                          minMarkerDist, "bp_", winName, "_",
                          propName, "_TelCenMatrix.txt"))
# Load TelCenMatrix
TelCenDF <- read.table(paste0("./cMMb_iwgsc_refseqv1.0_minInterMarkerDist",
                              minMarkerDist, "bp_", winName, "_",
                              propName, "_TelCenMatrix.txt"))
TelCenProfile <- as.vector(rowMeans(TelCenDF, na.rm = T))
TelCenSD <- as.vector(apply(X = TelCenDF, MARGIN = 1, FUN = sd))


# Function to plot telomere to centromere (Tel-Cen)
# profile of log2(markChIPA/markControlA) 
TelCenPlot <- function(xplot,
                       profile,
                       proportions,
                       proportionsName,
                       profileColour,
                       Ylabel,
                       Ylim) {
  plot(xplot, profile, type = "l", lwd = 3, col = profileColour,
       ylim = Ylim,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "")
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(1, seq(10, proportions, by = 10)),
       labels = c(expression(italic("TEL")),
                  seq(10, proportions-10, by = 10),
                  expression(italic("CEN"))))
  mtext(side = 1, line = 2.1, cex = 1,
        text = paste0("Scaled windows (", proportionsName, ")"))
  axis(side = 2, at = pretty(c(profile)), cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.1, cex = 1, text = Ylabel, col = profileColour)
  abline(h = 0, lwd = 1.5, lty = 1)
  box(lwd = 1.5)
#  legend(legendLoc,
#         legend = legendLabs,
#         col = c(profileColour),
#         text.col = c(profileColour),
#         text.font = c(1),
#         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

pdf(paste0("./cMMb_iwgsc_refseqv1.0_minInterMarkerDist",
           minMarkerDist, "bp_", winName, "_",
           propName, "_TelCenProfile.pdf"),
    height = 3.5, width = 7)
par(mfrow = c(1, 1))
par(mar = c(3.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
TelCenPlot(xplot = 1:length(TelCenProfile),
           profile = TelCenProfile,
           proportions = prop,
           proportionsName = propName,
           profileColour = profileColour,
           Ylabel = bquote("cM/Mb"),
           Ylim = c(min(TelCenProfile, na.rm = T),
                    max(TelCenProfile, na.rm = T)))
dev.off()
