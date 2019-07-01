#!/applications/R/R-3.5.0/bin/Rscript

# Divide genes into quantiles based on mean winName-scaled
# recombination rate (cM/Mb) from promoters to terminators.
# Plot average coverage profiles for a given dataset
# around these genes separated into cM/Mb quantiles.

# Usage:
#./cMMb_geneQuantiles_avgProfileRibbon.R 100kb 1 4 ASY1_CS ASY_CS_Rep1_ChIP input H3_input_SRR6350669 both purple4 genes Genes 20bp 20 2kb '2 kb' 2000 'euchromatin' '0.38,0.96'

#winName <- "100kb"
#minMarkerDist <- 1
#quantiles <- 4
#markChIP <- "ASY1_CS"
#libNameChIP <- "ASY1_CS_Rep1_ChIP"
#markControl <- "input"
#libNameControl <- "H3_input_SRR6350669"
#align <- "both"
#colour <- "purple4"
#featureName <- "genes"
#featureNamePlot <- "Genes"
#binName <- "20bp"
#binSize <- 20
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#upstream <- 2000
#downstream <- 2000
#region <- "euchromatin"
## top left
##legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
##                                        split = ",")))
## bottom left
##legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
##                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))


args <- commandArgs(trailingOnly = T)
winName <- args[1]
minMarkerDist <- as.numeric(args[2])
quantiles <- as.numeric(args[3])
markChIP <- args[4]
libNameChIP <- args[5]
markControl <- args[6]
libNameControl <- args[7]
align <- args[8]
colour <- args[9]
featureName <- args[10]
featureNamePlot <- args[11]
binName <- args[12]
binSize <- as.numeric(args[13])
flankName <- args[14]
flankNamePlot <- args[15]
upstream <- as.numeric(args[16])
downstream <- as.numeric(args[16])
region <- args[17]
legendPos <- as.numeric(unlist(strsplit(args[18],
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
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
} else {
  stop("region is not euchromatin, heterochromatin or genomewide")
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
} else if(region == "genomewide") {
  maskGR <- GRanges()
} else {
  stop("region is not euchromatin, heterochromatin or genomewide")
}

# Load table of representative genes and convert into GRanges
genes <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA.gff3",
                    colClasses = c(NA,
                                   rep("NULL", 2),
                                   rep(NA, 2),
                                   "NULL", NA, "NULL", NA))
colnames(genes) <- c("chr", "start", "end", "strand", "geneID")
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = genes$strand,
                   geneID = genes$geneID)
genesGR <- genesGR[seqnames(genesGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = maskGR,
                                   subject = genesGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
genesGR <- genesGR[-subjectHits(mask_genes_overlap)]
# Extend gene boundaries to include promoters and terminators for calculation of
# winName-scaled recombination rate
genesGR <- GRanges(seqnames = seqnames(genesGR),
                   ranges = IRanges(start = start(genesGR)-1000,
                                    end = end(genesGR)+1000),
                   strand = strand(genesGR),
                   geneID = genesGR$geneID)
print(genesGR)

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

# Obtain winName-scaled cMMb values for each gene between promoter and terminator
# Where features overlap more than one winName window, calculate mean cMMb
gene_cMMb_overlaps <- findOverlaps(query = genesGR,
                                   subject = cMMbGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
gene_cMMb_overlapsList <- lapply(seq_along(genesGR), function(x) {
  subjectHits(gene_cMMb_overlaps)[queryHits(gene_cMMb_overlaps) == x]
})
## OR
#gene_cMMb_overlapsList <- getOverlaps(coordinates = genesGR,
#                                      segments = cMMbGR,
#                                      overlapType = "overlapping",
#                                      whichOverlaps = TRUE,
#                                      ignoreStrand = TRUE)
gene_cMMb <- sapply(gene_cMMb_overlapsList,
                    function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))
genesGR <- GRanges(genesGR,
                   geneID = genesGR$geneID,
                   cMMb = gene_cMMb)
genesGR <- genesGR[genesGR$cMMb != "NaN"]

# Separate genes into genomes
genomeNames <- c("A", "B", "D", "")
genomeNamesPlot <- c("A", "B", "D", "All")
genesGR_genomeList <- lapply(seq_along(genomeNames), function(z) {
  genesGR[grep(genomeNames[z],
               seqnames(genesGR))]
})

# Order genes by decreasing cM/Mb
# and divide into quantiles for average coverage profiling
sapply(seq_along(genesGR_genomeList), function(z) {
  genesDF <- data.frame(genesGR_genomeList[[z]])
  genesDF_ordered <- genesDF[order(genesDF$cMMb,
                                   decreasing = T,
                                   na.last = NA),]
  quantilesN <- round(dim(genesDF_ordered)[1]/quantiles)
  quantilesCum <- cumsum(c(1,
                           rep(quantilesN,
                               times = quantiles)))
  if(quantilesCum[length(quantilesCum)] < dim(genesDF_ordered)[2]) {
    quantilesCum <- c(quantilesCum, dim(genesDF_ordered)[1])
  }
  quantilesStats <- data.frame()
  for(j in 1:length(quantilesCum)) {
    print(j)
    if(j == 1) {
      quantilejFeatures <- genesDF_ordered[(quantilesCum[j]):(quantilesCum[j+1]-1),]
      print(paste0("condition 1: quantile ", j))
      print(dim(quantilejFeatures))
    }
    if(j > 1 & j < quantiles) {
      quantilejFeatures <- genesDF_ordered[(quantilesCum[j]):(quantilesCum[j+1]-1),]
      print(paste0("condition 2: quantile ", j))
      print(dim(quantilejFeatures))
    }
    if(j == quantiles) {
      quantilejFeatures <- genesDF_ordered[(quantilesCum[j]):(quantilesCum[length(quantilesCum)]),]
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
                              featureName, "_in_", genomeNames[z], "genome_", region,
                              "_quantile", j, "_of_", quantiles, "quantiles_", 
                              "ordered_by_", winName, "Scaled_cMMb_",
                              "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    }
  }
  write.table(quantilesStats,
              file = paste0(regionDir,
                            "summary_", featureName, "_in_", genomeNames[z], "genome_", region, "_",
                            quantiles, "quantiles_", 
                            "ordered_by_", winName, "Scaled_cMMb_",
                            "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
})

# Load quantiles and order by geneID
genesDF_genomeList_quantiles <- lapply(seq_along(genesGR_genomeList),
  function(z) {
    lapply(seq_along(1:quantiles), function(j) {
      tmp <- read.table(paste0(regionDir,
                               featureName, "_in_", genomeNames[z], "genome_", region,
                               "_quantile", j, "_of_", quantiles, "quantiles_", 
                               "ordered_by_", winName, "Scaled_cMMb_",
                               "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                        header = T)
      #tmp[order(tmp$cMMb,
      #          decreasing = T,
      #          na.last = NA),]
      tmp[order(tmp$geneID,
                decreasing = F,
                na.last = NA),]
    })
})

# Load feature coverage matrices
## ChIP
if(libNameChIP %in% c("H3K4me3_ChIP_SRR6350668",
                      "H3K27me3_ChIP_SRR6350666",
                      "H3K36me3_ChIP_SRR6350670",
                      "H3K9ac_ChIP_SRR6350667",
                      "CENH3_ChIP_SRR1686799")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                       markChIP, "/snakemake_ChIPseq/mapped/geneProfiles/matrices/")
} else {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/",
                       markChIP, "/snakemake_ChIPseq/mapped/geneProfiles/matrices/")
}
covMatChIP <- read.table(paste0(covDirChIP,
                                libNameChIP,
                                "_MappedOn_wheat_v1.0_lowXM_",
                                align, "_sort_norm_",
                                featureName, "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         skip = 3)

## Control
if(libNameControl == "H3_input_SRR6350669") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                          "input/snakemake_ChIPseq/mapped/geneProfiles/matrices/")
} else if(libNameControl == "MNase_Rep1") {
  covDirControl <- paste0("/home/ajt200/analysis/wheat/",
                          "MNase/snakemake_ChIPseq/mapped/geneProfiles/matrices/")
} else {
  if(!(libNameControl %in% c("H3_input_SRR6350669", "MNase_Rep1"))) {
    stop("libNameControl is neither H3_input_SRR6350669 nor MNase_Rep1")
  }
}
covMatControl <- read.table(paste0(covDirControl,
                                   libNameControl,
                                   "_MappedOn_wheat_v1.0_lowXM_",
                                   align, "_sort_norm_",
                                   featureName, "_matrix_bin", binName,
                                   "_flank", flankName, ".tab"),
                            skip = 3)

## log2(ChIP/control) matrices
log2covMat <- log2((covMatChIP+1)/(covMatControl+1))

# Add geneIDs
log2covMat <- data.frame(geneID = genes$geneID,
                         log2covMat)
# Subdivide coverage matrix into above-defined quantiles
log2covMat_genomeList_quantiles <- lapply(seq_along(genesDF_genomeList_quantiles),
  function(z) {
    lapply(seq_along(1:quantiles), function(j) {
      tmp1 <- log2covMat[log2covMat$geneID %in% genesDF_genomeList_quantiles[[z]][[j]]$geneID,]
      tmp1[order(tmp1$geneID,
                 decreasing = F,
                 na.last = NA),]
    })
})
# Order by decreasing cM/Mb [THIS DOESN'T WORK - WORK OUT WHY]
log2covMat_genomeList_quantiles <- lapply(seq_along(log2covMat_genomeList_quantiles),
  function(z) {
    lapply(seq_along(1:quantiles), function(j) {
      log2covMat_genomeList_quantiles[[z]][[j]][order(genesDF_genomeList_quantiles[[z]][[j]]$cMMb,
                                                      decreasing = T,
                                                      na.last = NA),]
    })
})

## Check if ordering of geneIDs is consistent after
## ordering by decreasing cM/Mb
#sapply(seq_along(genesDF_genomeList_quantiles),
#  function(z) {
#    lapply(seq_along(1:quantiles), function(j) {
#      print(identical(as.character(genesDF_genomeList_quantiles[[z]][[j]]$geneID), as.character(log2covMat_genomeList_quantiles[[z]][[j]]$geneID)))
#    })
#})



## feature
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFfeature_list <- mclapply(seq_along(log2covMat_genomeList_quantiles), function(z) {
  lapply(seq_along(1:quantiles), function(j) {    
    data.frame(window = colnames(as.matrix(log2covMat_genomeList_quantiles[[z]][[j]][,-1])),
               t(as.matrix(log2covMat_genomeList_quantiles[[z]][[j]][,-1])))
    })
}, mc.cores = length(log2covMat_genomeList_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list  <- mclapply(seq_along(wideDFfeature_list), function(z) {
  lapply(seq_along(1:quantiles), function(j) {
    gather(data  = wideDFfeature_list[[z]][[j]],
           key   = feature,
           value = coverage,
           -window)
  })
}, mc.cores = length(wideDFfeature_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(z in seq_along(tidyDFfeature_list)) {
  for(j in seq_along(1:quantiles)) {
    tidyDFfeature_list[[z]][[j]]$window <- factor(tidyDFfeature_list[[z]][[j]]$window,
                                                  levels = as.character(wideDFfeature_list[[z]][[j]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list <- mclapply(seq_along(tidyDFfeature_list), function(z) {
  lapply(seq_along(1:quantiles), function(j) {
    data.frame(window = as.character(wideDFfeature_list[[z]][[j]]$window),
               n      = tapply(X     = tidyDFfeature_list[[z]][[j]]$coverage,
                               INDEX = tidyDFfeature_list[[z]][[j]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list[[z]][[j]]$coverage,
                               INDEX = tidyDFfeature_list[[z]][[j]]$window,
                               FUN   = mean),
               sd     = tapply(X     = tidyDFfeature_list[[z]][[j]]$coverage,
                               INDEX = tidyDFfeature_list[[z]][[j]]$window,
                               FUN   = sd))
  })
}, mc.cores = length(tidyDFfeature_list))

for(z in seq_along(summaryDFfeature_list)) {
  for(j in seq_along(1:quantiles)) {
    summaryDFfeature_list[[z]][[j]]$window <- factor(summaryDFfeature_list[[z]][[j]]$window,
                                                     levels = as.character(wideDFfeature_list[[z]][[j]]$window))
    summaryDFfeature_list[[z]][[j]]$winNo <- factor(1:dim(summaryDFfeature_list[[z]][[j]])[1])
    summaryDFfeature_list[[z]][[j]]$sem <- summaryDFfeature_list[[z]][[j]]$sd/sqrt(summaryDFfeature_list[[z]][[j]]$n-1)
    summaryDFfeature_list[[z]][[j]]$CI_lower <- summaryDFfeature_list[[z]][[j]]$mean -
      qt(0.975, df = summaryDFfeature_list[[z]][[j]]$n-1)*summaryDFfeature_list[[z]][[j]]$sem
    summaryDFfeature_list[[z]][[j]]$CI_upper <- summaryDFfeature_list[[z]][[j]]$mean +
      qt(0.975, df = summaryDFfeature_list[[z]][[j]]$n-1)*summaryDFfeature_list[[z]][[j]]$sem
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
for(z in seq_along(summaryDFfeature_list)) {
  names(summaryDFfeature_list[[z]]) <- quantileNames
}

# For each genome, convert list of lists summaryDFfeature_list[[z]]
# into a single data.frame of quantiles for plotting
summaryDFfeature_genome <- lapply(seq_along(summaryDFfeature_list), function(z) {
  bind_rows(summaryDFfeature_list[[z]], .id = "quantile")
})
for(z in seq_along(summaryDFfeature_genome)) {
  summaryDFfeature_genome[[z]]$quantile <- factor(summaryDFfeature_genome[[z]]$quantile,
                                                  levels = names(summaryDFfeature_list[[z]]))
}


## Random gene quantiles
# Define function to randomly select n rows from
# a GRanges object (features)
selectRandomFeatures <- function(features, n) {
  return(features[sample(x = dim(features)[1],
                         size = n,
                         replace = FALSE),])
}

# Define seed so that random selections are reproducible
set.seed(93750174)

# Divide genes into random sets of equal number,
# with the same number of genes per chromosome as
# above-defined cM/Mb-ordered gene quantiles 
sapply(seq_along(genesDF_genomeList_quantiles), function(z) {
  genesDF <- data.frame(genesGR_genomeList[[z]])
  randomPCsStats <- data.frame()
  for(j in 1:quantiles) {
    print(j)
    randomPCjFeatures <- data.frame()
    for(i in 1:length(chrs)) {
      randomPCjFeaturesChr <- selectRandomFeatures(features = genesDF[genesDF$seqnames
                                                                      == chrs[i],],
                                                   n = dim(genesDF_genomeList_quantiles[[z]][[j]][genesDF_genomeList_quantiles[[z]][[j]]$seqnames
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
                              featureName, "_in_", genomeNames[z], "genome_", region,
                              "_quantile", j, "_of_", quantiles, "quantiles_randomPerChrom_",
                              winName, "Scaled_cMMb_",
                              "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
  write.table(randomPCsStats,
              file = paste0(regionDir,
                            "summary_", featureName, "_in_", genomeNames[z], "genome_", region, "_",
                            quantiles, "quantiles_randomPerChrom_",
                            winName, "Scaled_cMMb_",
                            "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
})

# Load randomPCs and order by geneID
genesDF_genomeList_randomPCs <- lapply(seq_along(genesGR_genomeList),
  function(z) {
    lapply(seq_along(1:quantiles), function(j) {
      tmp <- read.table(paste0(regionDir,
                               featureName, "_in_", genomeNames[z], "genome_", region,
                               "_quantile", j, "_of_", quantiles, "quantiles_randomPerChrom_",
                               winName, "Scaled_cMMb_",
                               "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                        header = T)
      #tmp[order(tmp$cMMb,
      #          decreasing = T,
      #          na.last = NA),]
      tmp[order(tmp$geneID,
                decreasing = F,
                na.last = NA),]
    })
})

# Subdivide coverage matrix into above-defined randomPCs
log2covMat_genomeList_randomPCs <- lapply(seq_along(genesDF_genomeList_randomPCs),
  function(z) {
    lapply(seq_along(1:quantiles), function(j) {
      tmp1 <- log2covMat[log2covMat$geneID %in% genesDF_genomeList_randomPCs[[z]][[j]]$geneID,]
      tmp1[order(tmp1$geneID,
                 decreasing = F,
                 na.last = NA),]
    })
})


## randomF
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFrandomF_list <- mclapply(seq_along(log2covMat_genomeList_randomPCs), function(z) {
  lapply(seq_along(1:quantiles), function(j) {
    data.frame(window = colnames(as.matrix(log2covMat_genomeList_randomPCs[[z]][[j]][,-1])),
               t(as.matrix(log2covMat_genomeList_randomPCs[[z]][[j]][,-1])))
    })
}, mc.cores = length(log2covMat_genomeList_randomPCs))

# Convert into tidy data.frame (long format)
tidyDFrandomF_list  <- mclapply(seq_along(wideDFrandomF_list), function(z) {
  lapply(seq_along(1:quantiles), function(j) {
    gather(data  = wideDFrandomF_list[[z]][[j]],
           key   = randomF,
           value = coverage,
           -window)
  })
}, mc.cores = length(wideDFrandomF_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(z in seq_along(tidyDFrandomF_list)) {
  for(j in seq_along(1:quantiles)) {
    tidyDFrandomF_list[[z]][[j]]$window <- factor(tidyDFrandomF_list[[z]][[j]]$window,
                                                  levels = as.character(wideDFrandomF_list[[z]][[j]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (randomFs) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFrandomF_list <- mclapply(seq_along(tidyDFrandomF_list), function(z) {
  lapply(seq_along(1:quantiles), function(j) {
    data.frame(window = as.character(wideDFrandomF_list[[z]][[j]]$window),
               n      = tapply(X     = tidyDFrandomF_list[[z]][[j]]$coverage,
                               INDEX = tidyDFrandomF_list[[z]][[j]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFrandomF_list[[z]][[j]]$coverage,
                               INDEX = tidyDFrandomF_list[[z]][[j]]$window,
                               FUN   = mean),
               sd     = tapply(X     = tidyDFrandomF_list[[z]][[j]]$coverage,
                               INDEX = tidyDFrandomF_list[[z]][[j]]$window,
                               FUN   = sd))
  })
}, mc.cores = length(tidyDFrandomF_list))

for(z in seq_along(summaryDFrandomF_list)) {
  for(j in seq_along(1:quantiles)) {
    summaryDFrandomF_list[[z]][[j]]$window <- factor(summaryDFrandomF_list[[z]][[j]]$window,
                                                     levels = as.character(wideDFrandomF_list[[z]][[j]]$window))
    summaryDFrandomF_list[[z]][[j]]$winNo <- factor(1:dim(summaryDFrandomF_list[[z]][[j]])[1])
    summaryDFrandomF_list[[z]][[j]]$sem <- summaryDFrandomF_list[[z]][[j]]$sd/sqrt(summaryDFrandomF_list[[z]][[j]]$n-1)
    summaryDFrandomF_list[[z]][[j]]$CI_lower <- summaryDFrandomF_list[[z]][[j]]$mean -
      qt(0.975, df = summaryDFrandomF_list[[z]][[j]]$n-1)*summaryDFrandomF_list[[z]][[j]]$sem
    summaryDFrandomF_list[[z]][[j]]$CI_upper <- summaryDFrandomF_list[[z]][[j]]$mean +
      qt(0.975, df = summaryDFrandomF_list[[z]][[j]]$n-1)*summaryDFrandomF_list[[z]][[j]]$sem
  }
}

randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(z in seq_along(summaryDFrandomF_list)) {
  names(summaryDFrandomF_list[[z]]) <- randomPCNames
}

# For each genome, convert list of lists summaryDFrandomF_list[[z]]
# into a single data.frame of quantiles for plotting
summaryDFrandomF_genome <- lapply(seq_along(summaryDFrandomF_list), function(z) {
  bind_rows(summaryDFrandomF_list[[z]], .id = "quantile")
})
for(z in seq_along(summaryDFrandomF_genome)) {
  summaryDFrandomF_genome[[z]]$quantile <- factor(summaryDFrandomF_genome[[z]]$quantile,
                                                  levels = names(summaryDFrandomF_list[[z]]))
}


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
ymin <- min(sapply(seq_along(summaryDFfeature_genome), function(z) {
  min(c(summaryDFfeature_genome[[z]]$CI_lower),
      c(summaryDFrandomF_genome[[z]]$CI_lower))
}))
ymax <- max(sapply(seq_along(summaryDFfeature_genome), function(z) {
  max(c(summaryDFfeature_genome[[z]]$CI_upper),
      c(summaryDFrandomF_genome[[z]]$CI_upper))
}))

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
for(z in seq_along(summaryDFfeature_genome)) {
  ggObjGA <- NULL
  ggObj1 <- NULL
  ggObj1 <- ggplot(data = summaryDFfeature_genome[[z]],
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = quantile),
                  ) +
    geom_line(data = summaryDFfeature_genome[[z]],
              mapping = aes(colour = quantile),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFfeature_genome[[z]],
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
                                (dim(summaryDFfeature_list[[z]][[1]])[1])-(downstream/binSize),
                                dim(summaryDFfeature_list[[z]][[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                featureStartLab,
                                featureEndLab,
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFfeature_list[[z]][[1]])[1])-(downstream/binSize)),
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
    ggtitle(bquote(.(genomeNamesPlot[z]) * "-genome" ~ .(region) ~ 
                   .(featureName) ~ "(" * italic("n") ~ "=" ~
                   # INEXACT - NEED TO CHANGE
                   .(prettyNum(summaryDFfeature_genome[[z]]$n[1]*quantiles,
                               big.mark = ",", trim = T)) *
                   ")"))
  ggObj2 <- NULL
  ggObj2 <- ggplot(data = summaryDFrandomF_genome[[z]],
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = quantile),
                  ) +
    geom_line(data = summaryDFrandomF_genome[[z]],
              mapping = aes(colour = quantile),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFrandomF_genome[[z]],
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
                                (dim(summaryDFrandomF_list[[z]][[1]])[1])-(downstream/binSize),
                                dim(summaryDFrandomF_list[[z]][[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                featureStartLab,
                                featureEndLab,
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFrandomF_list[[z]][[1]])[1])-(downstream/binSize)),
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
    ggtitle(bquote(.(genomeNamesPlot[z]) * "-genome" ~ .(region) ~ 
                   .(featureName) ~ "(" * italic("n") ~ "=" ~
                   # INEXACT - NEED TO CHANGE
                   .(prettyNum(summaryDFrandomF_genome[[z]]$n[1]*quantiles,
                               big.mark = ",", trim = T)) *
                   ")"))
  ggObjGA <- grid.arrange(ggObj1, ggObj2, nrow = 1, ncol = 2)
  ggsave(paste0(plotDir,
                "log2_", libNameChIP, "_", libNameControl,
                "_around_", genomeNames[z], "genome_", region, "_",
                featureName, "_", quantiles, "quantiles_ordered_by_",
                "ordered_by_", winName, "Scaled_cMMb_",
                "minInterMarkerDist", as.character(minMarkerDist), "bp.pdf"),
         plot = ggObjGA,
         height = 8, width = 16)
}
