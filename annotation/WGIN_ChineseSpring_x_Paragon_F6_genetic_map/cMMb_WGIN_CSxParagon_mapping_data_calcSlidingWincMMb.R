#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./cMMb_WGIN_CSxParagon_mapping_data_calcSlidingWincMMb.R 10Mb 10000000 1Mb 1000000 1

#winName <- "10Mb"
#winSize <- 10000000
#stepName <- "1Mb"
#stepSize <- 1000000
#minMarkerDist <- 1

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
stepName <- args[3]
stepSize <- as.numeric(args[4])
minMarkerDist <- as.numeric(args[5])

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
  
  # Keep second marker in a pair of markers with an inter-marker distance >= minMarkerDist bp
  mapChrMMD <-  mapChr[1,]
  for(x in 2:dim(mapChr)[1]) {
    mapChrRow <- mapChr[x,][mapChr$physicalPosition[x]-mapChr$physicalPosition[x-1] >= minMarkerDist,]
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

print(paste0("Number of markers with physical distances >= ",
             minMarkerDist, " bp"))
print(dim(map_cMMb)[1])

windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  # Define sliding windows of width winSize nt,
  # with a step of stepSize nt
  ## Note: the commented-out code would create windows of winSize nt only,
  ## whereas the active code creates windows decreasing from winSize nt to stepSize nt
  ## at the right-hand end of each chromosome ( from chrLens[i]-winSize to chrLens[i] ),
  winStarts <- seq(from = 1,
                   to = chrLens[i],
#                   to = chrLens[i]-winSize,
                   by = stepSize)
#  if(chrLens[i] - winStarts[length(winStarts)] >= winSize) {
#    winStarts <- c(winStarts,
#                   winStarts[length(winStarts)]+stepSize)
#  }
  winEnds <- seq(from = winStarts[1]+winSize-1,
                 to = chrLens[i],
                 by = stepSize)
  winEnds <- c(winEnds,
               rep(chrLens[i], times = length(winStarts)-length(winEnds)))
  stopifnot(length(winStarts) == length(winEnds))

  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = IRanges(start = winStarts,
                                           end = winEnds),
                          strand = "*")
  print(chrWindowsGR)
  try( if( end(chrWindowsGR)[length(chrWindowsGR)] != chrLens[i] )
       {
         stop("End coordinate of last window != chromosome length")
       }
     )
  windowsGR <- append(windowsGR, chrWindowsGR)
}

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
                          as.character(minMarkerDist), "bp_", winName, "_step", stepName, ".txt"),
            sep = "\t", quote = F, row.names = F, col.names = T)
