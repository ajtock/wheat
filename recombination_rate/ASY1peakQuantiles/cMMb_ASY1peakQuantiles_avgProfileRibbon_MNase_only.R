#!/applications/R/R-3.5.0/bin/Rscript

# Divide peaks into quantiles based on mean winName-scaled
# recombination rate (cM/Mb) from promoters to terminators.
# Plot average coverage profiles for a given dataset
# around these peaks separated into cM/Mb quantiles.

# Usage:
#./cMMb_ASY1peakQuantiles_avgProfileRibbon_MNase_only.R 100kb 1 4 MNase MNase_Rep1 both purple4 ASY1_CS_peaks 'ASY1 peaks' 20bp 20 2kb '2 kb' 2000 'A' 'euchromatin' '0.02,0.94'

#winName <- "100kb"
#minMarkerDist <- 1
#quantiles <- 4
#markChIP <- "MNase"
#libNameChIP <- "MNase_Rep1"
#align <- "both"
#colour <- "purple4"
#featureName <- "ASY1_CS_peaks"
#featureNamePlot <- "ASY1 peaks"
#binName <- "20bp"
#binSize <- 20
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#upstream <- 2000
#downstream <- 2000
#genomeName <- "A"
#region <- "euchromatin"
## top left
##legendPos <- as.numeric(unlist(strsplit("0.02,0.94",
##                                        split = ",")))
## bottom left
##legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
##                                        split = ",")))
## top centre
##legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
##                                        split = ",")))


args <- commandArgs(trailingOnly = T)
winName <- args[1]
minMarkerDist <- as.numeric(args[2])
quantiles <- as.numeric(args[3])
markChIP <- args[4]
libNameChIP <- args[5]
align <- args[6]
colour <- args[7]
featureName <- args[8]
featureNamePlot <- args[9]
binName <- args[10]
binSize <- as.numeric(args[11])
flankName <- args[12]
flankNamePlot <- args[13]
upstream <- as.numeric(args[14])
downstream <- as.numeric(args[14])
genomeName <- args[15]
region <- args[16]
legendPos <- as.numeric(unlist(strsplit(args[17],
                                        split = ",")))

library(GenomicRanges)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

regionDir <- paste0(region, "/")
plotDir <- paste0(regionDir, "plots/")
system(paste0("[ -d ", regionDir, " ] || mkdir ", regionDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr2AMiddleInterval_chr4ALeftmostInterval_chr4BRightmostInterval_chr5ARightmostInterval_chr7BRightTwoIntervals.txt")[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
} else if(region == "pericentromeres") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
} else if(region == "pericentromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
} else if(region == "genomewide") {
  maskGR <- GRanges()
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Load peaks in BED format 
# Note addition of 1 to 0-based BED start coordinates
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           "ASY1_CS_Rep1_ChIP",
                           "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth_in_",
                           genomeName, "genome_", region, ".bed"),
                    header = F)
colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand")
peaksGR <- GRanges(seqnames = peaks$chr,
                   ranges = IRanges(start = peaks$start+1,
                                    end = peaks$end),
                   strand = peaks$strand,
                   peakID = peaks$name)
peaksGR <- peaksGR[seqnames(peaksGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
if(length(mask_peaks_overlap) > 0) {
 peaksGR <- peaksGR[-subjectHits(mask_peaks_overlap)]
}
# Extend peak boundaries to include promoters and terminators for calculation of
# winName-scaled recombination rate
peaksGR <- GRanges(seqnames = seqnames(peaksGR),
                   ranges = IRanges(start = start(peaksGR)-1000,
                                    end = end(peaksGR)+1000),
                   strand = strand(peaksGR),
                   peakID = peaksGR$peakID)
print(peaksGR)

# Convert windowed recombination rate into GRanges
cMMb <- read.table(paste0(
                   "/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/",
                   "cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist",
                   as.character(minMarkerDist), "bp_", winName, ".txt"))
cMMbGR <- GRanges(seqnames = cMMb$chr,
                  ranges = IRanges(start = cMMb$windowStart,
                                   end = cMMb$windowEnd),
                  strand = "*",
                  cMMb = cMMb$cMMb)

# Obtain winName-scaled cMMb values for each peak between promoter and terminator
# Where features overlap more than one winName window, calculate mean cMMb
peak_cMMb_overlaps <- findOverlaps(query = peaksGR,
                                   subject = cMMbGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peak_cMMb_overlapsList <- lapply(seq_along(peaksGR), function(x) {
  subjectHits(peak_cMMb_overlaps)[queryHits(peak_cMMb_overlaps) == x]
})
## OR
#peak_cMMb_overlapsList <- getOverlaps(coordinates = peaksGR,
#                                      segments = cMMbGR,
#                                      overlapType = "overlapping",
#                                      whichOverlaps = TRUE,
#                                      ignoreStrand = TRUE)
peak_cMMb <- sapply(peak_cMMb_overlapsList,
                    function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
peaksGR <- GRanges(peaksGR,
                   peakID = peaksGR$peakID,
                   cMMb = peak_cMMb)
peaksGR <- peaksGR[peaksGR$cMMb != "NaN"]

# Order peaks by decreasing cM/Mb
# and divide into quantiles for average coverage profiling
peaksDF <- data.frame(peaksGR)
peaksDF_ordered <- peaksDF[order(peaksDF$cMMb,
                                 decreasing = T,
                                 na.last = NA),]
quantilesN <- round(dim(peaksDF_ordered)[1]/quantiles)
quantilesCum <- cumsum(c(1,
                         rep(quantilesN,
                             times = quantiles)))
if(quantilesCum[length(quantilesCum)] < dim(peaksDF_ordered)[2]) {
  quantilesCum <- c(quantilesCum, dim(peaksDF_ordered)[1])
}
quantilesStats <- data.frame()
for(j in 1:length(quantilesCum)) {
  print(j)
  if(j == 1) {
    quantilejFeatures <- peaksDF_ordered[(quantilesCum[j]):(quantilesCum[j+1]-1),]
    print(paste0("condition 1: quantile ", j))
    print(dim(quantilejFeatures))
  }
  if(j > 1 & j < quantiles) {
    quantilejFeatures <- peaksDF_ordered[(quantilesCum[j]):(quantilesCum[j+1]-1),]
    print(paste0("condition 2: quantile ", j))
    print(dim(quantilejFeatures))
  }
  if(j == quantiles) {
    quantilejFeatures <- peaksDF_ordered[(quantilesCum[j]):(quantilesCum[length(quantilesCum)]),]
    print(paste0("condition 3: quantile ", j))
    print(dim(quantilejFeatures))
  }
  if(j <= quantiles) {
    stats <- data.frame(quantile = as.integer(j),
                        n = as.integer(dim(quantilejFeatures)[1]),
                        mean_width = as.integer(round(mean(quantilejFeatures$width, na.rm = T))),
                        total_width = as.integer(sum(quantilejFeatures$width, na.rm = T)),
                        mean_cMMb = as.numeric(mean(quantilejFeatures$cMMb, na.rm = T)))
  quantilesStats <- rbind(quantilesStats, stats)
  write.table(quantilejFeatures,
              file = paste0(regionDir,
                            featureName, "_in_", genomeName, "genome_", region,
                            "_quantile", j, "_of_", quantiles, "quantiles_", 
                            "ordered_by_", winName, "Scaled_cMMb_",
                            "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  }
}
write.table(quantilesStats,
            file = paste0(regionDir,
                          "summary_", featureName, "_in_", genomeName, "genome_", region, "_",
                          quantiles, "quantiles_", 
                          "ordered_by_", winName, "Scaled_cMMb_",
                          "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Load quantiles and order by peakID
peaksDF_quantiles <- lapply(seq_along(1:quantiles), function(j) {
  tmp <- read.table(paste0(regionDir,
                           featureName, "_in_", genomeName, "genome_", region,
                           "_quantile", j, "_of_", quantiles, "quantiles_", 
                           "ordered_by_", winName, "Scaled_cMMb_",
                           "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                    header = T)
  #tmp[order(tmp$cMMb,
  #          decreasing = T,
  #          na.last = NA),]
  tmp[order(tmp$peakID,
            decreasing = F,
            na.last = NA),]
})

# Load feature coverage matrices
## ChIP
if(libNameChIP %in% c("H3K4me3_ChIP_SRR6350668",
                      "H3K27me3_ChIP_SRR6350666",
                      "H3K36me3_ChIP_SRR6350670",
                      "H3K9ac_ChIP_SRR6350667",
                      "CENH3_ChIP_SRR1686799")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                       markChIP, "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
} else {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/",
                       markChIP, "/snakemake_ChIPseq/mapped/ASY1peakProfiles/matrices/")
}
covMatChIP <- read.table(paste0(covDirChIP,
                                libNameChIP,
                                "_MappedOn_wheat_v1.0_lowXM_",
                                align, "_sort_norm_",
                                featureName, "_in_", genomeName, "genome_", region,
                                "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         skip = 3)

# Add peakIDs
covMatChIP <- data.frame(peakID = peaks$name,
                         covMatChIP)
# Subdivide coverage matrix into above-defined quantiles
covMatChIP_quantiles <- lapply(seq_along(1:quantiles), function(j) {
   tmp1 <- covMatChIP[covMatChIP$peakID %in% peaksDF_quantiles[[j]]$peakID,]
   tmp1[order(tmp1$peakID,
              decreasing = F,
              na.last = NA),]
 })
# Order by decreasing cM/Mb [THIS DOESN'T WORK - WORK OUT WHY]
covMatChIP_quantiles <- lapply(seq_along(1:quantiles), function(j) {
  covMatChIP_quantiles[[j]][order(peaksDF_quantiles[[j]]$cMMb,
                                  decreasing = T,
                                  na.last = NA),]
})

## Check if ordering of peakIDs is consistent after
## ordering by decreasing cM/Mb
#lapply(seq_along(1:quantiles), function(j) {
#  print(identical(as.character(peaksDF_quantiles[[j]]$peakID),
#                  as.character(covMatChIP_quantiles[[j]]$peakID)))
#})



## feature
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFfeature_list <- lapply(seq_along(1:quantiles), function(j) {    
data.frame(window = colnames(as.matrix(covMatChIP_quantiles[[j]][,-1])),
           t(as.matrix(covMatChIP_quantiles[[j]][,-1])))
})

# Convert into tidy data.frame (long format)
tidyDFfeature_list  <- lapply(seq_along(1:quantiles), function(j) {
  gather(data  = wideDFfeature_list[[j]],
         key   = feature,
         value = coverage,
         -window)
})

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(j in seq_along(1:quantiles)) {
  tidyDFfeature_list[[j]]$window <- factor(tidyDFfeature_list[[j]]$window,
                                           levels = as.character(wideDFfeature_list[[j]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list <- lapply(seq_along(1:quantiles), function(j) {
  data.frame(window = as.character(wideDFfeature_list[[j]]$window),
             n      = tapply(X     = tidyDFfeature_list[[j]]$coverage,
                             INDEX = tidyDFfeature_list[[j]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFfeature_list[[j]]$coverage,
                             INDEX = tidyDFfeature_list[[j]]$window,
                             FUN   = mean,
                             na.rm = T),
             sd     = tapply(X     = tidyDFfeature_list[[j]]$coverage,
                             INDEX = tidyDFfeature_list[[j]]$window,
                             FUN   = sd,
                             na.rm = T))
})

for(j in seq_along(1:quantiles)) {
  summaryDFfeature_list[[j]]$window <- factor(summaryDFfeature_list[[j]]$window,
                                              levels = as.character(wideDFfeature_list[[j]]$window))
  summaryDFfeature_list[[j]]$winNo <- factor(1:dim(summaryDFfeature_list[[j]])[1])
  summaryDFfeature_list[[j]]$sem <- summaryDFfeature_list[[j]]$sd/sqrt(summaryDFfeature_list[[j]]$n-1)
  summaryDFfeature_list[[j]]$CI_lower <- summaryDFfeature_list[[j]]$mean -
    qt(0.975, df = summaryDFfeature_list[[j]]$n-1)*summaryDFfeature_list[[j]]$sem
  summaryDFfeature_list[[j]]$CI_upper <- summaryDFfeature_list[[j]]$mean +
    qt(0.975, df = summaryDFfeature_list[[j]]$n-1)*summaryDFfeature_list[[j]]$sem
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
names(summaryDFfeature_list) <- quantileNames

# For each genome, convert list of lists summaryDFfeature_list
# into a single data.frame of quantiles for plotting
summaryDFfeature <- bind_rows(summaryDFfeature_list, .id = "quantile")
summaryDFfeature$quantile <- factor(summaryDFfeature$quantile,
                                    levels = names(summaryDFfeature_list))


## Random peak quantiles
# Define function to randomly select n rows from
# a GRanges object (features)
selectRandomFeatures <- function(features, n) {
  return(features[sample(x = dim(features)[1],
                         size = n,
                         replace = FALSE),])
}

# Define seed so that random selections are reproducible
set.seed(93750174)

# Divide peaks into random sets of equal number,
# with the same number of peaks per chromosome as
# above-defined cM/Mb-ordered peak quantiles 
peaksDF <- data.frame(peaksGR)
randomPCsStats <- data.frame()
for(j in 1:quantiles) {
  print(j)
  randomPCjFeatures <- data.frame()
  for(i in 1:length(chrs)) {
    randomPCjFeaturesChr <- selectRandomFeatures(features = peaksDF[peaksDF$seqnames
                                                                    == chrs[i],],
                                                 n = dim(peaksDF_quantiles[[j]][peaksDF_quantiles[[j]]$seqnames
                                                                                == chrs[i],])[1])
    randomPCjFeatures <- rbind(randomPCjFeatures, randomPCjFeaturesChr)
  }
  print(dim(randomPCjFeatures))
  stats <- data.frame(quantile = as.integer(j),
                      n = as.integer(dim(randomPCjFeatures)[1]),
                      mean_width = as.integer(round(mean(randomPCjFeatures$width, na.rm = T))),
                      total_width = as.integer(sum(randomPCjFeatures$width, na.rm = T)),
                      mean_cMMb = as.numeric(mean(randomPCjFeatures$cMMb, na.rm = T)))
  randomPCsStats <- rbind(randomPCsStats, stats)
  write.table(randomPCjFeatures,
              file = paste0(regionDir,
                            featureName, "_in_", genomeName, "genome_", region,
                            "_quantile", j, "_of_", quantiles, "quantiles_randomPerChrom_",
                            winName, "Scaled_cMMb_",
                            "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}
write.table(randomPCsStats,
            file = paste0(regionDir,
                          "summary_", featureName, "_in_", genomeName, "genome_", region, "_",
                          quantiles, "quantiles_randomPerChrom_",
                          winName, "Scaled_cMMb_",
                          "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Load randomPCs and order by peakID
peaksDF_randomPCs <- lapply(seq_along(1:quantiles), function(j) {
  tmp <- read.table(paste0(regionDir,
                           featureName, "_in_", genomeName, "genome_", region,
                           "_quantile", j, "_of_", quantiles, "quantiles_randomPerChrom_",
                           winName, "Scaled_cMMb_",
                           "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                    header = T)
  #tmp[order(tmp$cMMb,
  #          decreasing = T,
  #          na.last = NA),]
  tmp[order(tmp$peakID,
            decreasing = F,
            na.last = NA),]
})

# Subdivide coverage matrix into above-defined randomPCs
covMatChIP_randomPCs <- lapply(seq_along(1:quantiles), function(j) {
  tmp1 <- covMatChIP[covMatChIP$peakID %in% peaksDF_randomPCs[[j]]$peakID,]
  tmp1[order(tmp1$peakID,
             decreasing = F,
             na.last = NA),]
})


## randomF
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFrandomF_list <- lapply(seq_along(1:quantiles), function(j) {
data.frame(window = colnames(as.matrix(covMatChIP_randomPCs[[j]][,-1])),
           t(as.matrix(covMatChIP_randomPCs[[j]][,-1])))
})

# Convert into tidy data.frame (long format)
tidyDFrandomF_list  <- lapply(seq_along(1:quantiles), function(j) {
  gather(data  = wideDFrandomF_list[[j]],
         key   = randomF,
         value = coverage,
         -window)
})

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(j in seq_along(1:quantiles)) {
  tidyDFrandomF_list[[j]]$window <- factor(tidyDFrandomF_list[[j]]$window,
                                           levels = as.character(wideDFrandomF_list[[j]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (randomFs) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFrandomF_list <- lapply(seq_along(1:quantiles), function(j) {
  data.frame(window = as.character(wideDFrandomF_list[[j]]$window),
             n      = tapply(X     = tidyDFrandomF_list[[j]]$coverage,
                             INDEX = tidyDFrandomF_list[[j]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFrandomF_list[[j]]$coverage,
                             INDEX = tidyDFrandomF_list[[j]]$window,
                             FUN   = mean,
                             na.rm = T),
             sd     = tapply(X     = tidyDFrandomF_list[[j]]$coverage,
                             INDEX = tidyDFrandomF_list[[j]]$window,
                             FUN   = sd,
                             na.rm = T))
})

for(j in seq_along(1:quantiles)) {
  summaryDFrandomF_list[[j]]$window <- factor(summaryDFrandomF_list[[j]]$window,
                                              levels = as.character(wideDFrandomF_list[[j]]$window))
  summaryDFrandomF_list[[j]]$winNo <- factor(1:dim(summaryDFrandomF_list[[j]])[1])
  summaryDFrandomF_list[[j]]$sem <- summaryDFrandomF_list[[j]]$sd/sqrt(summaryDFrandomF_list[[j]]$n-1)
  summaryDFrandomF_list[[j]]$CI_lower <- summaryDFrandomF_list[[j]]$mean -
    qt(0.975, df = summaryDFrandomF_list[[j]]$n-1)*summaryDFrandomF_list[[j]]$sem
  summaryDFrandomF_list[[j]]$CI_upper <- summaryDFrandomF_list[[j]]$mean +
    qt(0.975, df = summaryDFrandomF_list[[j]]$n-1)*summaryDFrandomF_list[[j]]$sem
}

randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
names(summaryDFrandomF_list) <- randomPCNames

# For each genome, convert list of lists summaryDFrandomF_list
# into a single data.frame of quantiles for plotting
summaryDFrandomF <- bind_rows(summaryDFrandomF_list, .id = "quantile")
summaryDFrandomF$quantile <- factor(summaryDFrandomF$quantile,
                                    levels = names(summaryDFrandomF_list))

# Specify x-axis tick labels
if(featureName == "genes" | featureName == "NLRs") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Define y-axis limits
#ymin <- min(c(summaryDFfeature$mean-summaryDFfeature$sem,
#              summaryDFranLoc$mean-summaryDFranLoc$sem))
#ymax <- max(c(summaryDFfeature$mean+summaryDFfeature$sem,
#              summaryDFranLoc$mean+summaryDFranLoc$sem))
ymin <- min(c(summaryDFfeature$CI_lower),
            c(summaryDFrandomF$CI_lower))
ymax <- max(c(summaryDFfeature$CI_upper),
            c(summaryDFrandomF$CI_upper))

# Function for formatting y-axis labels
# with a given number of decimals
fmt_decimals <- function(decimals) {
  function(x) format(x, nsmall = decimals, scientific = FALSE)
}

colours <- colorRampPalette(c("red", "blue"))(4)

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = colours[x], fontsize = 18)))
})
legendLabs_randomF <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = colours[x], fontsize = 18)))
})

# Plot average coverage profiles with 95% CI ribbon
## feature
ggObjGA <- NULL
ggObj1 <- NULL
ggObj1 <- ggplot(data = summaryDFfeature,
                 mapping = aes(x = winNo,
                               y = mean,
                               group = quantile),
                ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = colours) +
  geom_ribbon(data = summaryDFfeature,
              #mapping = aes(ymin = mean-sem,
              #              ymax = mean+sem,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(ymin, ymax),
                     labels = function(x) sprintf("%7.4f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_list[[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * "(" * .(markChIP) * "/control)")) +
  annotation_custom(legendLabs_feature[[1]]) +
  annotation_custom(legendLabs_feature[[2]]) +
  annotation_custom(legendLabs_feature[[3]]) +
  annotation_custom(legendLabs_feature[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        #legend.text = element_text(size = 10),
        #legend.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(colour = "transparent",
        #                          fill = "transparent"),
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle(bquote(.(genomeName) * "-genome" ~ .(region) ~ 
                 .(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 # INEXACT - NEED TO CHANGE
                 .(prettyNum(summaryDFfeature$n[1]*quantiles,
                             big.mark = ",", trim = T)) *
                 ")"))
ggObj2 <- NULL
ggObj2 <- ggplot(data = summaryDFrandomF,
                 mapping = aes(x = winNo,
                               y = mean,
                               group = quantile),
                ) +
  geom_line(data = summaryDFrandomF,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = colours) +
  geom_ribbon(data = summaryDFrandomF,
              #mapping = aes(ymin = mean-sem,
              #              ymax = mean+sem,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(ymin, ymax),
                     labels = function(x) sprintf("%7.4f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFrandomF_list[[1]])[1])-(downstream/binSize),
                              dim(summaryDFrandomF_list[[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFrandomF_list[[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * "(" * .(markChIP) * "/control)")) +
  annotation_custom(legendLabs_randomF[[1]]) +
  annotation_custom(legendLabs_randomF[[2]]) +
  annotation_custom(legendLabs_randomF[[3]]) +
  annotation_custom(legendLabs_randomF[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        #legend.text = element_text(size = 10),
        #legend.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(colour = "transparent",
        #                          fill = "transparent"),
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle(bquote(.(genomeName) * "-genome" ~ .(region) ~ 
                 .(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 # INEXACT - NEED TO CHANGE
                 .(prettyNum(summaryDFrandomF$n[1]*quantiles,
                             big.mark = ",", trim = T)) *
                 ")"))
ggObjGA <- grid.arrange(ggObj1, ggObj2, nrow = 1, ncol = 2)
ggsave(paste0(plotDir,
              libNameChIP, "_",
              "around_", genomeName, "genome_", region, "_",
              featureName, "_", quantiles, "quantiles_ordered_by_",
              "ordered_by_", winName, "Scaled_cMMb_",
              "minInterMarkerDist", as.character(minMarkerDist), "bp.pdf"),
       plot = ggObjGA,
       height = 7, width = 16)
