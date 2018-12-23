#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./cMMb_chrProfiles.R 1Mb 1000000

winName <- "1Mb"
winSize <- 1000000

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])

library(parallel)
library(doParallel)
library(GenomicRanges)
library(segmentSeq)

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_recombination_rate_analysis/"
outDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/"

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,3])

# Read in genetic map data
map <- read.table(paste0(inDir, "iwgsc_refseqv1.0_mapping_data.txt"),
                  header = T)

# Create numeric vectors of physical and genetic distances between markers
physicalDelta <- NULL
geneticDelta <- NULL
for(i in 1:length(chrs)) {
  mapChr <- map[map$chromosome == chrs[i],]
  
  physicalDeltaChr <- as.numeric(sapply(seq(from = 0,
                                            to = dim(mapChr)[1]-1),
    function(x) {
      mapChr$physicalPosition[x+1]-mapChr$physicalPosition[x]
    }))
  physicalDelta <- c(physicalDelta, physicalDeltaChr)

  geneticDeltaChr <- as.numeric(sapply(seq(from = 0,
                                           to = dim(mapChr)[1]-1),
    function(x) {
      mapChr$geneticPosition[x+1]-mapChr$geneticPosition[x]
    }))
  geneticDelta <- c(geneticDelta, geneticDeltaChr)
}

map_cMMb <- data.frame(map,
                       physicalDelta = as.numeric(physicalDelta),
                       geneticDelta = as.numeric(geneticDelta),
                       cMMb = as.numeric(log2((geneticDelta/(physicalDelta/1000000))+1)))
## Remove row at beginning of each chromosome with NA delta values (21 rows)
#map_cMMb <- map_cMMb[!is.na(map_cMMb$physicalDelta),]
# Note: there are three marker pairs with negative genetic distances

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
             window = as.integer(start(chrWindowsGR)),
             cMMb = win_cMMb)
}, mc.cores = length(chrs))

profile <- NULL
for(i in 1:length(chrProfiles)) {
  profile <- rbind(profile, chrProfiles[[i]])
}
write.table(profile,
            file = paste0(outDir,
                          "log2_cMMb_iwgsc_refseqv1.0_mapping_data_", winName, ".txt"))
