#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./cMMb_iwgsc_refseqv1.0_mapping_data_removeProblemMarkers_calcAdjacentWincMMb.R 1Mb 1000000 200 15

#winName <- "1Mb"
#winSize <- 1000000
#minMarkerDist <- 200
#N <- 15

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
minMarkerDist <- as.numeric(args[3])
N <- as.numeric(args[4])

library(parallel)
library(doParallel)
library(GenomicRanges)
library(segmentSeq)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)

# Read in genetic and physical map data
map <- read.table("iwgsc_refseqv1.0_mapping_data.txt",
                  header = T)
print("Number of pre-filtered markers:")
print(dim(map)[1])

# Create numeric vectors of physical and genetic distances between markers
mapMMD <- NULL
physicalDelta <- NULL
geneticDelta <- NULL
for(i in 1:length(chrs)) {
  mapChr <- map[map$chromosome == chrs[i],]
  
  # Keep second marker in a pair of markers with an inter-marker distance >= minMarkerDist bp
  # and with a genetic distance >= 0 cM
  mapChrMMD <-  mapChr[1,]
  for(x in 2:dim(mapChr)[1]) {
    mapChrRow <- mapChr[x,][mapChr$physicalPosition[x]-mapChr$physicalPosition[x-1] >= minMarkerDist &
                            mapChr$geneticPosition[x]-mapChr$geneticPosition[x-1] >= 0,]
    mapChrMMD <- rbind(mapChrMMD, mapChrRow)
  }
  mapMMD <- rbind(mapMMD, mapChrMMD)

  physicalDeltaChr <- as.numeric(sapply(seq(from = 1,
                                            to = dim(mapChrMMD)[1]),
    function(x) {
      mapChrMMD$physicalPosition[x]-mapChrMMD$physicalPosition[x-1]
    }
  ))
  physicalDelta <- c(physicalDelta, physicalDeltaChr)

  geneticDeltaChr <- as.numeric(sapply(seq(from = 1,
                                           to = dim(mapChrMMD)[1]),
    function(x) {
      mapChrMMD$geneticPosition[x]-mapChrMMD$geneticPosition[x-1]
    }
  ))
  geneticDelta <- c(geneticDelta, geneticDeltaChr)
}

map_cMMb <- data.frame(mapMMD,
                       physicalDelta    = as.numeric(physicalDelta),
                       geneticDelta     = as.numeric(geneticDelta),
                       cMMb             = as.numeric(geneticDelta/(physicalDelta/1000000)))

print(paste0("Number of markers with inter-marker genetic distances >= 0 and physical distances >= ",
             minMarkerDist, " bp"))
print(dim(map_cMMb)[1])

# Remove second marker for 557 marker pairs with negative genetic distances
map_cMMb <- map_cMMb[is.na(map_cMMb$geneticDelta) | map_cMMb$geneticDelta >= 0,]
print("Number of markers with positive inter-marker genetic distances")
print(dim(map_cMMb)[1])

## Remove second marker for pairs of markers with physical distances < minMarkerDist bp
#map_cMMb <- map_cMMb[is.na(map_cMMb$physicalDelta) | map_cMMb$physicalDelta >= minMarkerDist,]

windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  seqWindows <- seq(from = 1,
                    to = chrLens[i],
                    by = winSize)
  windowsIR <- IRanges(start = seqWindows,
                       width = winSize)
  windowsIR <- windowsIR[-length(windowsIR)]
  windowsIR <- append(windowsIR,
                      IRanges(start = seqWindows[length(seqWindows)],
                              end = chrLens[i]))
  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = windowsIR,
                          strand = "*")
  print(chrWindowsGR)
  try( if( end(chrWindowsGR)[length(chrWindowsGR)] != chrLens[i] )
       {
         stop("End coordinate of last window != chromosome length")
       }
     )
  windowsGR <- append(windowsGR, chrWindowsGR)
}

registerDoParallel(cores = length(chrs))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

chrProfiles <- mclapply(seq_along(chrs), function(i) {
  map_cMMbChr <- map_cMMb[map_cMMb$chromosome == chrs[i],]
  markerIntervalsGR <- GRanges()
  for(k in 1:(dim(map_cMMbChr)[1]-1)) {
    markerIntervalsGRrow <- GRanges(seqnames = chrs[i],
                                    ranges = IRanges(start = map_cMMbChr$physicalPosition[k],
                                                     end = map_cMMbChr$physicalPosition[k+1]),
                                    cMMb = map_cMMbChr$cMMb[k+1])
    markerIntervalsGR <- append(markerIntervalsGR, markerIntervalsGRrow)
  }
   
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  overlaps <- getOverlaps(chrWindowsGR,
                          markerIntervalsGR,
                          whichOverlaps = T,
                          ignoreStrand = T)
  win_cMMb <- sapply(overlaps,
                     function(x) mean(markerIntervalsGR$cMMb[x]))
  data.frame(chr = as.character(seqnames(chrWindowsGR)),
             windowStart = as.integer(start(chrWindowsGR)),
             windowEnd = as.integer(end(chrWindowsGR)),
             cMMb = win_cMMb)
}, mc.cores = length(chrs))

profile <- NULL
for(i in 1:length(chrProfiles)) {
  profile <- rbind(profile, chrProfiles[[i]])
}
write.table(profile,
            file = paste0("cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                          as.character(minMarkerDist), "bp_", winName, ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

cMMbProfile <- read.table(paste0("cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                                 as.character(minMarkerDist), "bp_", winName, ".txt"),
                          header = T)

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  cMMbProfile[cMMbProfile$chr == chrs[x],]
}, mc.cores = length(chrs))

filt_chrProfiles <- mclapply(seq_along(chrProfiles), function(x) {
  filt_chrProfile <- stats::filter(x = chrProfiles[[x]]$cMMb,
                                   filter = f,
                                   sides = 2)
  # Given missing cM/Mb data for some of the more distal windows,
  # need a different way of extending the leftmost and rightmost
  # non-NA values to the ends of each chromosome, replacing NAs where they are present
  leftFlank <- which(is.na(filt_chrProfile))[which(is.na(filt_chrProfile)) < N*2]
  rightFlank <- which(is.na(filt_chrProfile))[which(is.na(filt_chrProfile)) > N*2]
  filt_chrProfile[leftFlank] <- filt_chrProfile[leftFlank[length(leftFlank)]+1]
  filt_chrProfile[rightFlank] <- filt_chrProfile[rightFlank[1]-1]
#  filt_chrProfile[1:flank] <- filt_chrProfile[flank+1]
#  filt_chrProfile[(length(filt_chrProfile)-flank+1):length(filt_chrProfile)] <- filt_chrProfile[(length(filt_chrProfile)-flank)]
  data.frame(chr = as.character(chrProfiles[[x]]$chr),
             windowStart = as.integer(chrProfiles[[x]]$windowStart),
             windowEnd = as.integer(chrProfiles[[x]]$windowEnd),
             filt_cMMb = as.numeric(filt_chrProfile),
             stringsAsFactors = F)
}, mc.cores = length(chrProfiles))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_cMMbProfile <- do.call(rbind, filt_chrProfiles)

write.table(filt_cMMbProfile,
            file = paste0("cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                          as.character(minMarkerDist), "bp_", winName, "_smoothed.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
