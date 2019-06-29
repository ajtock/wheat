#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./coverage_per_scaled_win_TelCen_noZscore.R H3K27me3 H3K27me3_ChIP_SRR6350666 input H3_input_SRR6350669 both 1Mb 1000000 100 100ths navy

#markChIPA <- "H3K27me3"
#libNameChIPA <- "H3K27me3_ChIP_SRR6350666"
#markControlA <- "input"
#libNameControlA <- "H3_input_SRR6350669"
#align <- "both"
#winName <- "1Mb"
#winSize <- 1000000
#prop <- 100
#propName <- "100ths" 
#profileColour <- "navy"

args <- commandArgs(trailingOnly = T)
markChIPA <- args[1]
libNameChIPA <- args[2]
markControlA <- args[3]
libNameControlA <- args[4]
align <- args[5]
winName <- args[6]
winSize <- as.numeric(args[7])
prop <- as.numeric(args[8])
propName <- args[9]
profileColour <- args[10]

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

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
centromere <- centromereStart+((centromereEnd-centromereStart)/2)

## ChIPA profile
if(libNameChIPA %in% c("H3K4me3_ChIP_SRR6350668",
                       "H3K27me3_ChIP_SRR6350666",
                       "H3K36me3_ChIP_SRR6350670",
                       "H3K9ac_ChIP_SRR6350667",
                       "CENH3_ChIP_SRR1686799")) {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
} else {
  covDirChIPA <- paste0("/home/ajt200/analysis/wheat/",
                        markChIPA, "/snakemake_ChIPseq/mapped/",
                        align, "/bg/")
}

profileChIPA <- read.table(paste0(covDirChIPA, libNameChIPA, "_MappedOn_wheat_v1.0_lowXM_",
                                  align, "_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileChIPA_bigWins <- profileChIPA[profileChIPA$V3-profileChIPA$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileChIPA <- profileChIPA[profileChIPA$V3-profileChIPA$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end 
profileChIPA_bigWinsList <- mclapply(seq_along(1:dim(profileChIPA_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileChIPA_bigWins[x,]$V2,
                      to = profileChIPA_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileChIPA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPA_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileChIPA_bigWins[x,]$V3)),
               V4 = as.numeric(profileChIPA_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileChIPA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileChIPA_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileChIPA_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileChIPA_bigWinsDT <- rbindlist(profileChIPA_bigWinsList)
profileChIPA <- rbind.fill(profileChIPA, profileChIPA_bigWinsDT)
profileChIPA <- profileChIPA[order(profileChIPA$V1, profileChIPA$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileChIPA <- profileChIPA[profileChIPA$V1 == chrs[x],]
  if(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileChIPA[dim(chrProfileChIPA)[1],]$V4))
  }
}, mc.cores = detectCores())
profileChIPA_chrLenValsDT <- rbindlist(chrLenValsList)
profileChIPA <- rbind.fill(profileChIPA, profileChIPA_chrLenValsDT)
profileChIPA <- profileChIPA[order(profileChIPA$V1, profileChIPA$V2),]

profileChIPA <- data.frame(chr = as.character(profileChIPA$V1),
                           window = as.integer(profileChIPA$V2+1),
                           CPM = as.numeric(profileChIPA$V4),
                           stringsAsFactors = F)

## ControlA profile
if(libNameControlA %in% c("H3_input_SRR6350669")) {
  covDirControlA <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                            markControlA, "/snakemake_ChIPseq/mapped/",
                            align, "/bg/")
} else {
  covDirControlA <- paste0("/home/ajt200/analysis/wheat/",
                            markControlA, "/snakemake_ChIPseq/mapped/",
                            align, "/bg/")
}

profileControlA <- read.table(paste0(covDirControlA, libNameControlA, "_MappedOn_wheat_v1.0_lowXM_",
                                     align, "_sort_norm_binSize", winName, ".bedgraph"))
# Rows where the difference between end and start coordinates is > winSize
profileControlA_bigWins <- profileControlA[profileControlA$V3-profileControlA$V2 > winSize,]
# Rows where the difference between end and start coordinates is == winSize
profileControlA <- profileControlA[profileControlA$V3-profileControlA$V2 == winSize,]

# Create a list of big windows, each split into windows of winSize,
# or < winSize if at chromosome end 
profileControlA_bigWinsList <- mclapply(seq_along(1:dim(profileControlA_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = profileControlA_bigWins[x,]$V2,
                      to = profileControlA_bigWins[x,]$V3,
                      by = winSize)

  if(bigWinsSplit[length(bigWinsSplit)] < profileControlA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlA_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+winSize,
                                 profileControlA_bigWins[x,]$V3)),
               V4 = as.numeric(profileControlA_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == profileControlA_bigWins[x,]$V3) {
    data.frame(V1 = as.character(profileControlA_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+winSize),
               V4 = as.numeric(profileControlA_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

profileControlA_bigWinsDT <- rbindlist(profileControlA_bigWinsList)
profileControlA <- rbind.fill(profileControlA, profileControlA_bigWinsDT)
profileControlA <- profileControlA[order(profileControlA$V1, profileControlA$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileControlA <- profileControlA[profileControlA$V1 == chrs[x],]
  if(chrProfileControlA[dim(chrProfileControlA)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileControlA[dim(chrProfileControlA)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileControlA[dim(chrProfileControlA)[1],]$V4))
  }
}, mc.cores = detectCores())
profileControlA_chrLenValsDT <- rbindlist(chrLenValsList)
profileControlA <- rbind.fill(profileControlA, profileControlA_chrLenValsDT)
profileControlA <- profileControlA[order(profileControlA$V1, profileControlA$V2),]

profileControlA <- data.frame(chr = as.character(profileControlA$V1),
                              window = as.integer(profileControlA$V2+1),
                              CPM = as.numeric(profileControlA$V4),
                              stringsAsFactors = F)

# Calculate log2((ChIPA+1)/(ControlA+1)) coverage within each window
profileLog2 <- data.frame(chr = as.character(profileChIPA$chr),
                          window = as.numeric(profileChIPA$window),
                          log2CPM = as.numeric(log2((profileChIPA$CPM+1)/(profileControlA$CPM+1))),
                          stringsAsFactors = F)
profileLog2 <- data.frame(profileLog2,
                          ZscoreLog2CPM = (profileLog2$log2CPM-mean(profileLog2$log2CPM,
                                                                    na.rm = T))/sd(profileLog2$log2CPM,
                                                                                   na.rm = T)) 
 
chrProfiles <- mclapply(seq_along(chrs), function(x) {
  profileLog2[profileLog2$chr == chrs[x],]
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

  # Calculate mean log2(markChIPA/markControlA) in each
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
                        mean(chrProfiles[[x]]$log2CPM[y])
                      })
  scaledWinAvgCovR <- sapply(overlapsR, function(y) {
                        mean(chrProfiles[[x]]$log2CPM[y])
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
            file = paste0("./log2_",
                          libNameChIPA, "_",
                          libNameControlA, "_",
                          align, "_", winName, "_",
                          propName, "_TelCenMatrix_noZscore.txt"))

# Load TelCenMatrix
TelCenDF <- read.table(paste0("./log2_",
                              libNameChIPA, "_",
                              libNameControlA, "_",
                              align, "_", winName, "_",
                              propName, "_TelCenMatrix_noZscore.txt"))
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

pdf(paste0("./log2_",
           libNameChIPA, "_",
           libNameControlA, "_",
           align, "_", winName, "_",
           propName, "_TelCenProfile_noZscore.pdf"),
    height = 3.5, width = 7)
par(mfrow = c(1, 1))
par(mar = c(3.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
TelCenPlot(xplot = 1:length(TelCenProfile),
           profile = TelCenProfile,
           proportions = prop,
           proportionsName = propName,
           profileColour = profileColour,
           Ylabel = bquote("Log"[2]*"(ChIP/input)"),
           Ylim = c(min(TelCenProfile),
                    max(TelCenProfile)),
           legendLabs = markChIPA,
           legendLoc = "top")
dev.off()
