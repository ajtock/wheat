#!/applications/R/R-3.5.0/bin/Rscript

# Divide genes into quantiles based on mean winName-scaled
# recombination rate (cM/Mb) from promoters to terminators.
# Plot average coverage profiles for a given dataset
# around these genes separated into cM/Mb quantiles.

# Usage:
#./cMMb_geneQuantiles_avgProfileRibbon.R 100kb 1 4 ASY1_CS ASY_CS_Rep1_ChIP input H3_input_SRR6350669 both purple4 genes 20bp 20 2kb 2000 'euchromatin'

winName <- "100kb"
minMarkerDist <- 1
quantiles <- 4
markChIP <- "ASY1_CS"
libNameChIP <- "ASY1_CS_Rep1_ChIP"
markControl <- "input"
libNameControl <- "H3_input_SRR6350669"
align <- "both"
colour <- "purple4"
featureName <- "genes"
binName <- "20bp"
binSize <- 20
flankName <- "2kb"
flankSize <- 2000
region <- "euchromatin"

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
binName <- args[11]
binSize <- as.numeric(args[12])
flankName <- args[13]
flankSize <- as.numeric(args[14])
region <- args[15]

library(GenomicRanges)
library(parallel)

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
genomeNamesPlot <- c("A", "B", "D", "all")
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
  if(quantilesCum[length(quantilesCum)] < dim(genesDF_ordered)[1]) {
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
      stats <-  data.frame(quantile = as.integer(j),
                           n = as.integer(dim(quantilejFeatures)[1]),
                           mean_width = as.integer(round(mean(quantilejFeatures$width))),
                           total_width = as.integer(sum(quantilejFeatures$width)),
                           mean_cMMb = as.numeric(mean(quantilejFeatures$cMMb)))
    quantilesStats <- rbind(quantilesStats, stats)
    write.table(quantilejFeatures,
                file = paste0(regionDir,
                              "genes_in_", genomeNames[z], "genome_", region,
                              "_quantile", j, "_of_", quantiles, "quantiles_", 
                              "ordered_by_", winName, "Scaled_cMMb_",
                              "minInterMarkerDist", as.character(minMarkerDist), "bp.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    }
  }
  write.table(quantilesStats,
              file = paste0(regionDir,
                            "summary_genes_in_", genomeNames[z], "genome_", region, "_",
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
                               "genes_in_", genomeNames[z], "genome_", region,
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

#sapply(seq_along(genesDF_genomeList_quantiles),
#  function(z) {
#    lapply(seq_along(1:quantiles), function(j) {
#      print(identical(as.character(genesDF_genomeList_quantiles[[z]][[j]]$geneID), as.character(log2covMat_genomeList_quantiles[[z]][[j]]$geneID)))
#    })
#})

