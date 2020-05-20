#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./TE_superfams_per_scaled_win_TelCen.R 1Mb 1000000 100 100ths grey50

#winName <- "1Mb"
#winSize <- 1000000
#prop <- 100
#propName <- "100ths" 
#profileColour <- "grey50"

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
prop <- as.numeric(args[3])
propName <- args[4]
profileColour <- args[5]

library(parallel)
library(segmentSeq)
library(doParallel)

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

# TE superfam
inDirSuperfams <- "/home/ajt200/analysis/wheat/chromosomeProfiles/TEs/superfamilies/"
superfamName <- c(
                  "Mariner_DTT",
                  "Mutator_DTM",
                  "Harbinger_DTH",
                  "CACTA_DTC",
                  "hAT_DTA",
                  "Unclassified_with_TIRs_DTX",
                  "Helitrons_DHH",
                  "Unclassified_class_2_DXX",
                  "LINE_RIX",
                  "SINE_SIX",
                  "Copia_LTR_RLC",
                  "Gypsy_LTR_RLG",
                  "Unclassified_LTR_RLX",
                  "Unclassified_repeats_XXX"
)
superfamNamePlot <- c(
                      "Mariner",
                      "Mutator",
                      "Harbinger",
                      "CACTA",
                      "hAT",
                      "Unclassified with TIRs",
                      "Helitron",
                      "Unclassified class 2",
                      "LINE",
                      "SINE",
                      "Copia LTR",
                      "Gypsy LTR",
                      "Unclassified LTR",
                      "Unclassified repeats"
)

superfamProfile <- mclapply(seq_along(superfamName), function(h) {
  read.table(paste0(inDirSuperfams,
                    "TE_frequency_per_", winName,
                    "_superfamily_", superfamName[h],
                    "_unsmoothed.txt"),
             header = T)
}, mc.cores = length(superfamName))

superfamChrProfiles <- mclapply(seq_along(superfamProfile), function(h) {
  lapply(seq_along(chrs), function(i) {
    superfamProfile[[h]][superfamProfile[[h]]$chr == chrs[i],]
  })
})

registerDoParallel(cores = length(superfamName))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Calculate average TEs in proportionally scaled
# windows along left and right chromosome arms
foreach(h = 1:length(superfamName)) %dopar% {
  print(superfamName[h])
  TelCenMatrix <- NULL
  for(x in seq_along(chrs)) {
    print(chrs[x])
    # Create proportionally scaled TEL-CEN windows
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
    # corresponding to winName windowed TEs values
    TEsWindowStarts <- superfamChrProfiles[[h]][[x]]$window
    TEsWindowEnds <- c((superfamChrProfiles[[h]][[x]]$window[2:length(superfamChrProfiles[[h]][[x]]$window)])-1,
                        chrLens[x])
    TEsWindowGRanges <- GRanges(seqnames = chrs[x],
                                ranges = IRanges(start = TEsWindowStarts,
                                                 end = TEsWindowEnds),
                                strand = "*")
  
    # Calculate mean TEs in each
    # scaled window using TEs values in winName windows
    overlapsL <- getOverlaps(coordinates = windowGRangesL,
                             segments = TEsWindowGRanges,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE)
    overlapsR <- getOverlaps(coordinates = windowGRangesR,
                             segments = TEsWindowGRanges,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE)
    scaledWinAvgCovL <- sapply(overlapsL, function(y) {
                          mean(superfamChrProfiles[[h]][[x]][,dim(superfamChrProfiles[[h]][[x]])[2]][y], na.rm = T)
                        })
    scaledWinAvgCovR <- sapply(overlapsR, function(y) {
                          mean(superfamChrProfiles[[h]][[x]][,dim(superfamChrProfiles[[h]][[x]])[2]][y], na.rm = T)
                        })
    scaledWinAvgCovLR <- cbind(scaledWinAvgCovL,
                               scaledWinAvgCovR)
    TelCenMatrix <- cbind(TelCenMatrix, scaledWinAvgCovLR)
  }
  write.table(TelCenMatrix,
              file = paste0("TEs_",
                            winName, "_",
                            "superfamily_", superfamName[h], "_",
                            propName, "_TelCenMatrix.txt"))
  # Load TelCenMatrix
  TelCenDF <- read.table(paste0("TEs_",
                                winName, "_",
                                "superfamily_", superfamName[h], "_",
                                propName, "_TelCenMatrix.txt"))
  TelCenProfile <- as.vector(rowMeans(TelCenDF, na.rm = T))
  TelCenSD <- as.vector(apply(X = TelCenDF, MARGIN = 1, FUN = sd, na.rm = T))
  
  
  # Function to plot telomere to centromere (Tel-Cen) profile
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
  
  pdf(paste0(plotDir, "TEs_",
             winName, "_",
             "superfamily_", superfamName[h], "_",
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
             Ylabel = superfamNamePlot[h],
             Ylim = c(min(TelCenProfile, na.rm = T),
                      max(TelCenProfile, na.rm = T)))
  dev.off()
}
