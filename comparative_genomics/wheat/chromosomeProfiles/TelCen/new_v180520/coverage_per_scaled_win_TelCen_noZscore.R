#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./coverage_per_scaled_win_TelCen_noZscore.R log2ChIPcontrol log2_DMC1_Rep1_ChIP_H3_input_SRR6350669 both 1Mb 1000000 100 100ths green2

#inDir <- "log2ChIPcontrol"
#profileName <- "log2_DMC1_Rep1_ChIP_H3_input_SRR6350669"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#prop <- 100
#propName <- "100ths" 
#profileColour <- "green2"

args <- commandArgs(trailingOnly = T)
inDir <- args[1]
profileName <- args[2]
align <- args[3]
winName <- args[4]
winSize <- as.numeric(args[5])
prop <- as.numeric(args[6])
propName <- args[7]
profileColour <- args[8]

library(parallel)
library(data.table)
library(plyr)
library(GenomicRanges)
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

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                               "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,2])
centromereEnd <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                             "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
centromere <- (centromereStart+centromereEnd)/2

# Load pre-calculated windowed coverage values or feature counts
profiles <- read.table(paste0("/home/ajt200/analysis/wheat/chromosomeProfiles/",
                              inDir, "/", profileName, "_per_", winName, "_unsmoothed.txt"),
                       header = T)

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  profiles[profiles$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate average coverage in proportionally scaled
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
  # Check that last window start coordinate == centromere[x]-1
  # (i.e., most proximal coordinate)
  stopifnot(end(windowGRangesL[length(windowGRangesL)]) == centromere[x]-1)
  seqWindowStartsR <- as.integer(seq(from = centromere[x],
                                     to = chrLens[x],
                                     by = ((chrLens[x]-centromere[x])+1)/prop))
  # Check that last window start coordinate is
  # < chrLens[x]
  stopifnot(seqWindowStartsR[length(seqWindowStartsR)] < chrLens[x])

  # Right arm
  seqWindowEndsR <- c(seqWindowStartsR[2:length(seqWindowStartsR)]-1,
                      chrLens[x])
  windowGRangesR <- rev(GRanges(seqnames = chrs[x],
                                 ranges = IRanges(start = seqWindowStartsR,
                                                  end = seqWindowEndsR),
                                 strand = "*"))
  print(windowGRangesR)
  # Check that first window start coordinate == chrLens[x]
  # (i.e., most distal coordinate)
  stopifnot(end(windowGRangesR[1]) == chrLens[x])

  # Create GRanges object of winName windows
  # corresponding to winName windowed coverage values
  covWindowStarts <- chrProfiles[[x]]$window
  covWindowEnds <- c((chrProfiles[[x]]$window[2:length(chrProfiles[[x]]$window)])-1,
                      chrLens[x])
  covWindowGRanges <- GRanges(seqnames = chrs[x],
                              ranges = IRanges(start = covWindowStarts,
                                               end = covWindowEnds),
                              strand = "*")

  # Calculate mean value in each
  # scaled window using coverage values in winName windows
  overlapsL <- getOverlaps(coordinates = windowGRangesL,
                           segments = covWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  overlapsR <- getOverlaps(coordinates = windowGRangesR,
                           segments = covWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  scaledWinAvgCovL <- sapply(overlapsL, function(y) {
                        mean(chrProfiles[[x]][,dim(chrProfiles[[x]])[2]][y], na.rm = T)
                      })
  scaledWinAvgCovR <- sapply(overlapsR, function(y) {
                        mean(chrProfiles[[x]][,dim(chrProfiles[[x]])[2]][y], na.rm = T)
                      })
  scaledWinAvgCovLR <- cbind(scaledWinAvgCovL,
                             scaledWinAvgCovR)
  TelCenMatrix <- cbind(TelCenMatrix, scaledWinAvgCovLR)
  #scaledWinAvgCovMeanLR <- sapply(seq_along(scaledWinAvgCovL), function(y) {
  #                           mean(c(scaledWinAvgCovL[y], scaledWinAvgCovR[y]))
  #                         })
  #scaledWinAvgCov <- scaledWinAvgCovL+scaledWinAvgCovR
  #TelCenProfile <- TelCenProfile+scaledWinAvgCov
}
write.table(TelCenMatrix,
            file = paste0(profileName, "_",
                          align, "_", winName, "_",
                          propName, "_TelCenMatrix.txt"))

# Load TelCenMatrix
TelCenDF <- read.table(paste0(profileName, "_",
                              align, "_", winName, "_",
                              propName, "_TelCenMatrix.txt"))
TelCenProfile <- as.vector(rowMeans(TelCenDF))
TelCenSD <- as.vector(apply(X = TelCenDF, MARGIN = 1, FUN = sd))


# Function to plot telomere to centromere (Tel-Cen) profile 
TelCenPlot <- function(xplot,
                       profile,
                       proportions,
                       proportionsName,
                       profileColour,
                       Ylabel,
                       Ylim,
                       legendLabs,
                       legendLoc) {
  plot(xplot, profile, type = "h", lwd = 4.9, col = profileColour,
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
  mtext(side = 2, line = 2.1, cex = 1, text = Ylabel, col = "black")
  abline(h = 0, lwd = 1.5, lty = 1)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(profileColour),
         text.col = c(profileColour),
         text.font = c(1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

pdf(paste0(plotDir, profileName,
           "_", align, "_", winName, "_",
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
           Ylabel = "Coverage",
#           Ylabel = bquote("Log"[2]*"(ChIP/input)"),
           Ylim = c(min(TelCenProfile),
                    max(TelCenProfile)),
           legendLabs = profileName,
           legendLoc = "top")
dev.off()
