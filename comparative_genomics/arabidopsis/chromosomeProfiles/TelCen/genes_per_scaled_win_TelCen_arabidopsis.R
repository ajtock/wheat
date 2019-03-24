#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./genes_per_scaled_win_TelCen_arabidopsis.R 10kb 10000 100 100ths forestgreen

winName <- "10kb"
winSize <- 10000
prop <- 100
propName <- "100ths" 
profileColour <- "forestgreen"

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
prop <- as.numeric(args[3])
propName <- args[4]
profileColour <- args[5]

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
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromere <- c(15086045, 3607929, 13587786, 3956021, 11725024)

## genes profile
inDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/genomewideZscore/"

genes <- read.table(paste0(inDir, "gene_frequency_genome_",
                           winName, ".txt"))

chrgeneprofiles <- mclapply(seq_along(chrs), function(x) {
  chrgenes <- genes[genes$chr == chrs[x],]
}, mc.cores = length(chrs))
chrProfiles <- mclapply(seq_along(chrgeneprofiles), function(x) {
  data.frame(chr    = as.character(chrgeneprofiles[[x]]$chr),
             window = seq(from = 1,
                          to = chrLens[x],
                          by = winSize),
             genes   = as.integer(chrgeneprofiles[[x]]$genes))
}, mc.cores = length(chrgeneprofiles))
             
# Calculate average genes in proportionally scaled
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
  # corresponding to winName windowed genes values
  genesWindowStarts <- chrProfiles[[x]]$window
  genesWindowEnds <- c((chrProfiles[[x]]$window[2:length(chrProfiles[[x]]$window)])-1,
                       chrLens[x])
  genesWindowGRanges <- GRanges(seqnames = chrs[x],
                               ranges = IRanges(start = genesWindowStarts,
                                                end = genesWindowEnds),
                               strand = "*")

  # Calculate mean genes in each
  # scaled window using genes values in winName windows
  overlapsL <- getOverlaps(coordinates = windowGRangesL,
                           segments = genesWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  overlapsR <- getOverlaps(coordinates = windowGRangesR,
                           segments = genesWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  scaledWinAvgCovL <- sapply(overlapsL, function(y) {
                        mean(chrProfiles[[x]]$genes[y], na.rm = T)
                      })
  scaledWinAvgCovR <- sapply(overlapsR, function(y) {
                        mean(chrProfiles[[x]]$genes[y], na.rm = T)
                      })
  scaledWinAvgCovLR <- cbind(scaledWinAvgCovL,
                             scaledWinAvgCovR)
  TelCenMatrix <- cbind(TelCenMatrix, scaledWinAvgCovLR)
}
write.table(TelCenMatrix,
            file = paste0("./genes_per_",
                          winName, "_",
                          propName, "_TelCenMatrix.txt"))
# Load TelCenMatrix
TelCenDF <- read.table(paste0("./genes_per_",
                              winName, "_",
                              propName, "_TelCenMatrix.txt"))
TelCenProfile <- as.vector(rowMeans(TelCenDF))
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

pdf(paste0("./genes_per_",
           winName, "_",
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
           Ylabel = bquote("Genes"),
           Ylim = c(min(TelCenProfile, na.rm = T),
                    max(TelCenProfile, na.rm = T)))
dev.off()
