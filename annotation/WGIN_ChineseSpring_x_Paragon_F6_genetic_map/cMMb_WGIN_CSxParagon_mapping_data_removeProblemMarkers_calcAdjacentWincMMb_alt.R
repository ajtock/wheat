#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./cMMb_WGIN_CSxParagon_mapping_data_removeProblemMarkers_calcAdjacentWincMMb_alt.R 1Mb 1000000 1 15

#winName <- "1Mb"
#winSize <- 1000000
#minMarkerDist <- 1
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
# "ChromRefSeq" and "SnpRefSeqPosition" give the physical marker positions
# Discard rows where "ChromRefSeq" (RefSeq v1.0 chromosome) does not match
# "ConsGenChrom" (genetic map chromosome) to exclude spuriously aligned markers
map <- read.csv("Brinton_Uauy_WGIN_CSxParagon_mapping_data.csv",
                header = T, stringsAsFactors = F)
print("Number of pre-filtered markers:")
print(dim(map)[1])
#[1] 35143

map <- map[,colnames(map) %in% c("Axiom35K_SNPId", "Axiom820K_Affycode",
                                 "ChromRefSeq", "SnpRefSeqPosition",
                                 "ConsGenChrom", "consensus_cM")]
map$ConsGenChrom <- sub("", "chr", map$ConsGenChrom)
map <- map[!is.na(map$ChromRefSeq) &
           !is.na(map$ConsGenChrom) &
           map$ChromRefSeq == map$ConsGenChrom,]
print("Number of post-filtered markers:")
print(dim(map)[1])
#[1] 9713

# Simplify table for calculation of cM/Mb
map <- map[,colnames(map) %in% c("Axiom35K_SNPId", "ChromRefSeq",
                                 "SnpRefSeqPosition", "consensus_cM")]
colnames(map) <- c("Axiom35K_SNPId", "chromosome", "physicalPosition", "geneticPosition")
map <- map[order(map$chromosome,
                 map$physicalPosition),]

# Create numeric vectors of physical and genetic distances between markers
mapMMD <- NULL
physicalDelta <- NULL
geneticDelta <- NULL
for(i in 1:length(chrs)) {
  mapChr <- map[map$chromosome == chrs[i],]
  mapChrMMD <-  mapChr[1,]
  # Initiate x at 2 to iterate over mapChr markers from
  # marker 2 to third-from-last marker ( dim(mapChr)[1]-1 )
  x <- 2
  while(x < (dim(mapChr)[1]-1)) {
    # Filter markers so that marker x is retained when the the following conditions are met:
    #                       the physical distance between each pair of markers is >= minMarkerDist AND
    #                   ( ( the genetic distance between between a marker [x] and the preceding marker [x-1] is >= 0 cM AND
    #                       the genetic distance between between the next marker [x+1] and the original marker [x] is >= 0 cM ) OR
    #                     ( the genetic distance between between a marker [x] and the preceding marker [x-1] is >= 0 cM AND
    #                       the genetic distance between between the next, next marker [x+2] and the original marker [x] is >= 0 cM ) )
    # The x+2 provision enables retention of marker x that has a cM position
    # greater than that of x+1 but less than or equal to that of x+2
    # I.e., rather than identifying and removing marker x as spurious,
    # this provision enables identification and removal of marker x+1 as spurious
    # where the cM position of x+1 is less than that of x AND greater than that of x+2 OR x+3
    mapChrRow <- mapChr[x,][mapChr$physicalPosition[x]-mapChr$physicalPosition[x-1] >= minMarkerDist &
                        ( ( mapChr$geneticPosition[x]-mapChr$geneticPosition[x-1] >= 0 &
                            mapChr$geneticPosition[x+1]-mapChr$geneticPosition[x] >= 0 ) |
                          ( mapChr$geneticPosition[x]-mapChr$geneticPosition[x-1] >= 0 &
                            mapChr$geneticPosition[x+2]-mapChr$geneticPosition[x] >= 0 ) ),]
    if(dim(mapChrRow)[1] != 0) {
      mapChrMMD <- rbind(mapChrMMD, mapChrRow)
      mapChr <- rbind(mapChrMMD, mapChr[(x+1):(dim(mapChr)[1]),])
      x <- x+1 
    } else {
      mapChr <- rbind(mapChrMMD, mapChr[(x+1):(dim(mapChr)[1]),])
    }
  }
  # Deal with final marker
  x <- dim(mapChr)[1]
  mapChrRow <- mapChr[x,][mapChr$physicalPosition[x]-mapChr$physicalPosition[x-1] >= minMarkerDist &
                          mapChr$geneticPosition[x]-mapChr$geneticPosition[x-1] >= 0,]
  if(dim(mapChrRow)[1] == 0) {
    mapChr <- mapChr[-x,]
  }
  mapMMD <- rbind(mapMMD, mapChr)

  physicalDeltaChr <- as.numeric(sapply(seq(from = 1,
                                            to = dim(mapChr)[1]),
    function(x) {
      mapChr$physicalPosition[x]-mapChr$physicalPosition[x-1]
    }
  ))
  physicalDelta <- c(physicalDelta, physicalDeltaChr)

  geneticDeltaChr <- as.numeric(sapply(seq(from = 1,
                                           to = dim(mapChr)[1]),
    function(x) {
      mapChr$geneticPosition[x]-mapChr$geneticPosition[x-1]
    }
  ))
  geneticDelta <- c(geneticDelta, geneticDeltaChr)
}

map_cMMb <- data.frame(mapMMD,
                       physicalDelta    = as.numeric(physicalDelta),
                       geneticDelta     = as.numeric(geneticDelta),
                       cMMb             = as.numeric(geneticDelta/(physicalDelta/1000000)))

print(paste0("Number of markers with physical distances >= ",
             minMarkerDist, " bp"))
print(dim(map_cMMb)[1])
#[1] 5946

# Remove second marker for 10 marker pairs with negative genetic distances
map_cMMb <- map_cMMb[is.na(map_cMMb$geneticDelta) | map_cMMb$geneticDelta >= 0,]
print("Number of markers with positive inter-marker genetic distances")
print(dim(map_cMMb)[1])
#[1] 5936

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
            file = paste0("cMMb_WGIN_CSxParagon_mapping_data_minInterMarkerDist",
                          as.character(minMarkerDist), "bp_", winName, "_alt.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)

# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

cMMbProfile <- read.table(paste0("cMMb_WGIN_CSxParagon_mapping_data_minInterMarkerDist",
                                 as.character(minMarkerDist), "bp_", winName, "_alt.txt"),
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
            file = paste0("cMMb_WGIN_CSxParagon_mapping_data_minInterMarkerDist",
                          as.character(minMarkerDist), "bp_", winName, "_alt_smoothed.txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
