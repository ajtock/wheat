#!/applications/R/R-3.5.0/bin/Rscript

#
# Divide features into quantiles based on population genetics statistics
# (e.g., Tajima's D, Pi normalised by locus with, or nucleotide diversity normalised by locus width).
# in a given feature region (e.g., TSS to TTS).
# Extract and save feature IDs for each quantile for further analyses
# (e.g., GO enrichment and average + 95% CI profile plotting).
# Calculate mean winName-scaled recombination rate (cM/Mb) from
# the promoter to the terminator of each feature.
#

# Usage:
# /applications/R/R-3.4.0/bin/Rscript group_genes_into_exomeSNP_quantiles.R 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' bodies 4 100kb 200

featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
                               split = ","))
region <- "bodies"
quantiles <- 4
winName <- "100kb"
minMarkerDist <- 200

args <- commandArgs(trailingOnly = T)
featureName <- unlist(strsplit(args[1],
                               split = ","))
region <- args[2]
quantiles <- as.numeric(args[3])
winName <- args[4]
minMarkerDist <-  as.numeric(args[5])

library(PopGenome)
library(WhopGenome)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

outDir <- paste0("quantiles_by_nucleotideDiversity_in_",
                 region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load features
features <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                    paste0(substring(featureName[x], first = 10, last = 16),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 18), ".gff3"),
             header = F)
})
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature data.frames
if(length(featureName) == 3) {
 features <- do.call(rbind, features)
} else {
 features <- features[[1]]
}
featureIDs <- sub(pattern = "\\.\\d+", replacement = "",
                  features$V9)

## Subset features to only those corresponding to NLRs
#features_NLRs <- features[features$V9 %in% NLRs$V9,]
## Get NLR IDs and their row indices in features
#IDs <- sub(pattern = "\\.\\d+", replacement = "",
#           features_NLRs$V9)
#ID_indices <- which(featureIDs %in% IDs)
#nonIDs <- featureIDs[!(featureIDs %in% IDs)]
## Function to randomly select feature IDs not present in IDs
#ran_nonIDs_select <- function(nonIDsChr, n) {
#  sample(x = nonIDsChr,
#         size = n,
#         replace = FALSE)
#}
## Apply ran_nonIDs_select() function on a per-chromosome basis
## and create growing vector of feature IDs called ran_nonIDs
#set.seed(9237452)
#ran_nonIDs <- NULL
#for(i in 1:length(levels(features$V1))) {
#  IDsChr <- IDs[grepl(paste0("TraesCS",
#                             sub("chr", "", levels(features$V1))[i]),
#                      IDs)]
#  nonIDsChr <- nonIDs[grepl(paste0("TraesCS",
#                                   sub("chr", "", levels(features$V1))[i]),
#                      nonIDs)]
#  ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
#                                     n = length(IDsChr))
#  ran_nonIDs <- c(ran_nonIDs, ran_nonIDsChr)
#}
#ran_nonID_indices <- which(featureIDs %in% ran_nonIDs)

# Load 811 exomes-derived SNP matrix from He et al. (2019) Nat. Genet.
if(!(file.exists("/home/ajt200/analysis/wheat/annotation/221118_download/He_Akhunov_2019_NatGenet_1000exomes_SNPs/all.GP08_mm75_het3_publication01142019.vcf.gz.tbi"))) {
  tabix_build("/home/ajt200/analysis/wheat/annotation/221118_download/He_Akhunov_2019_NatGenet_1000exomes_SNPs/all.GP08_mm75_het3_publication01142019.vcf.gz",
              sc = as.integer(1), bc = as.integer(2), ec = as.integer(2), meta = "#", lineskip = as.integer(0))
} else {
  print("*.vcf.gz.tbi index file exists")
}
vcf_handle <- vcf_open("/home/ajt200/analysis/wheat/annotation/221118_download/He_Akhunov_2019_NatGenet_1000exomes_SNPs/all.GP08_mm75_het3_publication01142019.vcf.gz")
#chrs <- unique(as.character(features$V1)
chrs <- paste0(rep("chr", 21), rep(1:7, 3),
               c(rep("A", 7), rep("B", 7), rep("D", 7)))
# This works with lapply, but fails with mclapply for some reason
genomeClass_list <- lapply(seq_along(chrs), function(x) {
  print(x)
  Whop_readVCF(v = vcf_handle,
               numcols = 1000,
               tid = chrs[x],
               frompos = 1,
               topos = 1000000000,
               samplenames = NA,
               gffpath = "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3",
               include.unknown = T)
})

# Extact variants located within features
genomeClassSplit_list <- mclapply(seq_along(genomeClass_list), function(x) {
  splitting.data(genomeClass_list[[x]],
                 positions = lapply(which(features$V1 == chrs[x]), function(x) {
                               features[x,]$V4:features[x,]$V5 }),
                 type = 2)
}, mc.cores = detectCores(), mc.preschedule = F)

# Get neutrality, diversity, and F_ST (fixation index) statistics
genomeClassSplit_list <- mclapply(seq_along(genomeClassSplit_list), function(x) {
  neutrality.stats(genomeClassSplit_list[[x]],
                   FAST = F, do.R2 = T)
}, mc.cores = detectCores(), mc.preschedule = F)
neutrality_stats_df_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
  data.frame(get.neutrality(genomeClassSplit_list[[x]],
                            theta = T)[[1]],
             stringsAsFactors = F)
})

genomeClassSplit_list <- mclapply(seq_along(genomeClassSplit_list), function(x) {
  diversity.stats(genomeClassSplit_list[[x]],
                  pi = T, keep.site.info = T)
}, mc.cores = detectCores(), mc.preschedule = F)
diversity_stats_df_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
  data.frame(get.diversity(genomeClassSplit_list[[x]],
                           between = F)[[1]],
             stringsAsFactors = F)
})

#genomeClassSplit_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
#  F_ST.stats(genomeClassSplit_list[[x]],
#             detail = T, mode = "nucleotide", FAST = F)
#})
## This generates only 0s or NaNs for NLRs
#F_ST_stats_df_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
#  data.frame(get.F_ST(genomeClassSplit_list[[x]],
#             mode = "nucleotide"),
#             stringsAsFactors = F)
#})

# Combine statistics into one dataframe
popgen_stats_df_list <- lapply(seq_along(neutrality_stats_df_list), function(x) {
  data.frame(chr = as.character(chrs[x]),
             start = as.integer(sub(pattern = " - \\d+", replacement = "",
                                    x = genomeClassSplit_list[[x]]@region.names)),
             end = as.integer(sub(pattern = "\\d+ - ", replacement = "",
                                  x = genomeClassSplit_list[[x]]@region.names)),
             width = as.integer( ( as.integer(sub(pattern = "\\d+ - ", replacement = "",
                                                  x = genomeClassSplit_list[[x]]@region.names)) -
                                   as.integer(sub(pattern = " - \\d+", replacement = "",
                                                  x = genomeClassSplit_list[[x]]@region.names)) ) + 1 ),
             ## OR
             #width = as.integer(genomeClassSplit_list[[x]]@n.sites),
             strand = as.character(features[features$V1 == chrs[x],]$V7),
             ID = as.character(features[features$V1 == chrs[x],]$V9),
             nuc.diversity.within = as.numeric(diversity_stats_df_list[[x]]$nuc.diversity.within),
             hap.diversity.within = as.numeric(diversity_stats_df_list[[x]]$hap.diversity.within),
             Pi = as.numeric(diversity_stats_df_list[[x]]$Pi),
             hap.F_ST.vs.all = as.numeric(diversity_stats_df_list[[x]]$hap.F_ST.vs.all),
             nuc.F_ST.vs.all = as.numeric(diversity_stats_df_list[[x]]$nuc.F_ST.vs.all),
             neutrality_stats_df_list[[x]],
             stringsAsFactors = F,
             row.names = as.character(1:dim(neutrality_stats_df_list[[x]])[1]))
})
popgen_stats <- do.call(rbind, popgen_stats_df_list)
# Sanity check to ensure popgen_stats feature chromosome IDs and start and end coordinates
# are identical to those in features data.frame
if(!identical(popgen_stats$chr, features$V1)) { 
  stop("popgen_stats chromosome IDs are not identical to those in features data.frame!")
}
if(!identical(popgen_stats$start, features$V4)) {
  stop("popgen_stats start coordinates are not identical to those in features data.frame!")
}
if(!identical(popgen_stats$end, features$V5)) {
  stop("popgen_stats end coordinates are not identical to those in features data.frame!")
}
    

# Load features 
features <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                    paste0(substring(featureName[x], first = 10, last = 16),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 18), ".gff3"),
             colClasses = c(NA,
                            rep("NULL", 2),
                            rep(NA, 2),
                            "NULL", NA, "NULL", NA),
             header = F)
})
if(length(featureName) == 3) {
  features <- do.call(rbind, features)
} else {
  features <- features[[1]]
}
colnames(features) <- c("chr", "start", "end", "strand", "featureID")
featuresGR <- GRanges(seqnames = features$chr,
                      ranges = IRanges(start = features$start,
                                       end = features$end),
                      strand = features$strand,
                      featureID = features$featureID)
# Extend feature boundaries to include promoters and terminators for calculation of
# winName-scaled recombination rate
featuresGR_ext <- GRanges(seqnames = seqnames(featuresGR),
                          ranges = IRanges(start = start(featuresGR)-1000,
                                           end = end(featuresGR)+1000),
                          strand = strand(featuresGR),
                          featureID = featuresGR$featureID)
print(featuresGR_ext)

# Load ranLocs 
ranLocs <- lapply(seq_along(featureName), function(y) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                    substring(featureName[y], first = 10), "_randomLoci.bed"),
             colClasses = c(rep(NA, 4), "NULL", NA),
             header = F)
})
if(length(featureName) == 3) {
  ranLocs <- do.call(rbind, ranLocs)
} else {
  ranLocs <- ranLocs[[1]]
}
colnames(ranLocs) <- c("chr", "start", "end", "ranLocID", "strand")
ranLocsGR <- GRanges(seqnames = ranLocs$chr,
                     ranges = IRanges(start = ranLocs$start+1,
                                      end = ranLocs$end),
                     strand = ranLocs$strand,
                     ranLocID = ranLocs$ranLocID)
# Extend feature boundaries to include promoters and terminators for calculation of
# winName-scaled recombination rate
ranLocsGR_ext <- GRanges(seqnames = seqnames(ranLocsGR),
                         ranges = IRanges(start = start(ranLocsGR)-1000,
                                          end = end(ranLocsGR)+1000),
                         strand = strand(ranLocsGR),
                         ranLocID = ranLocsGR$ranLocID)
print(ranLocsGR_ext)

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

# Obtain winName-scaled cMMb values for each feature between promoter and terminator
# Where features overlap more than one winName window, calculate mean cMMb
feature_cMMb_overlaps <- findOverlaps(query = featuresGR_ext,
                                      subject = cMMbGR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
feature_cMMb_overlapsList <- lapply(seq_along(featuresGR_ext), function(x) {
  subjectHits(feature_cMMb_overlaps)[queryHits(feature_cMMb_overlaps) == x]
})
## OR using segmentSeq::getOverlaps()
#feature_cMMb_overlapsList <- getOverlaps(coordinates = featuresGR_ext,
#                                         segments = cMMbGR,
#                                         overlapType = "overlapping",
#                                         whichOverlaps = TRUE,
#                                         ignoreStrand = TRUE)
feature_cMMb <- sapply(feature_cMMb_overlapsList,
                       function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))

# Get intron number per gene
introns <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                    "IWGSC_v1.1_HC_20170706_representative_introns_in_",
                    paste0(substring(featureName[x], first = 10, last = 16),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 18), ".gff3"),
             header = F,
             stringsAsFactors = F)
})
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature data.frames
if(length(featureName) == 3) {
  introns <- do.call(rbind, introns)
} else {
  introns <- introns[[1]]
}
colnames(introns) <- c("seqid", "source", "type", "start", "end",
                       "score", "strand", "phase", "ID", "Parent")
intronNoPerGene <- unlist(mclapply(seq_along(as.character(features$featureID)), function(x) {
  dim(introns[as.character(introns$Parent) == as.character(features$featureID)[x],])[1]
}, mc.cores = detectCores()))

# Get exon number per gene
exons <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                    "IWGSC_v1.1_HC_20170706_representative_exons_in_",
                    paste0(substring(featureName[x], first = 10, last = 16),
                           collapse = "_"), "_",
                    substring(featureName[1][1], first = 18), ".gff3"),
             header = F,
             stringsAsFactors = F)
})
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature data.frames
if(length(featureName) == 3) {
  exons <- do.call(rbind, exons)
} else {
  exons <- exons[[1]]
}
colnames(exons) <- c("seqid", "source", "type", "start", "end",
                     "score", "strand", "phase", "ID", "Parent")
exonNoPerGene <- unlist(mclapply(seq_along(as.character(features$featureID)), function(x) {
  dim(exons[as.character(exons$Parent) == as.character(features$featureID)[x],])[1]
}, mc.cores = detectCores()))

# Combine feature coordinates and their corresponding mean haplotype numbers,
# cM/Mb values, and intron and exon numbers
featuresGR <- GRanges(featuresGR,
                      featureID = featuresGR$featureID,
                      # Normalise by feature width (per 1 kb)
                      nucleotideDiversity = popgen_stats$nuc.diversity.within/(popgen_stats$width/1e3),
                      # The PopGenome developers advise that, where include.unknown = TRUE in the
                      # Whop_readVCF() function call, "diversity.stats: pi and haplotype diversity should not be used"
                      haplotypeDiversity = popgen_stats$hap.diversity.within,
                      # Normalise by feature width (per 1 kb)
                      Pi = popgen_stats$Pi/(popgen_stats$width/1e3),
                      TajimaD = popgen_stats$Tajima.D,
                      nSegregatingSites = popgen_stats$n.segregating.sites/(popgen_stats$width/1e3),
                      RozasR2 = popgen_stats$Rozas.R_2,
                      FuLiF = popgen_stats$Fu.Li.F,
                      FuLiD = popgen_stats$Fu.Li.D,
                      cMMb = feature_cMMb,
                      exons = exonNoPerGene,
                      introns = intronNoPerGene)
#featuresGR <- featuresGR[ID_indices]

# Sanity check to ensure featuresGR feauture IDs, chromosome IDs and start and end coordinates
# are identical to those in features data.frame
if(!identical(gsub("\\.\\d+", "", as.character(featuresGR$featureID)), featureIDs)) {
  stop("featuresGR feature IDs are not identical to those in featureIDs vector!")
}
if(!identical(as.character(seqnames(featuresGR)), as.character(features$chr))) { 
  stop("featuresGR chromosome IDs are not identical to those in features data.frame!")
}
if(!identical(start(featuresGR), features$start)) {
  stop("featuresGR start coordinates are not identical to those in features data.frame!")
}
if(!identical(end(featuresGR), features$end)) {
  stop("featuresGR end coordinates are not identical to those in features data.frame!")
}

# Obtain winName-scaled cMMb values for each ranLoc between promoter and terminator
# Where ranLocs overlap more than one winName window, calculate mean cMMb
ranLoc_cMMb_overlaps <- findOverlaps(query = ranLocsGR_ext,
                                     subject = cMMbGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
ranLoc_cMMb_overlapsList <- lapply(seq_along(ranLocsGR_ext), function(x) {
  subjectHits(ranLoc_cMMb_overlaps)[queryHits(ranLoc_cMMb_overlaps) == x]
})
## OR using segmentSeq::getOverlaps()
#ranLoc_cMMb_overlapsList <- getOverlaps(coordinates = ranLocsGR_ext,
#                                        segments = cMMbGR,
#                                        overlapType = "overlapping",
#                                        whichOverlaps = TRUE,
#                                        ignoreStrand = TRUE)
ranLoc_cMMb <- sapply(ranLoc_cMMb_overlapsList,
                      function(x) mean(cMMbGR$cMMb[x], na.rm = TRUE))

# Combine ranLoc coordinates and their corresponding mean haplotype numbers,
# cM/Mb values 
ranLocsGR <- GRanges(ranLocsGR,
                     ranLocID = ranLocsGR$ranLocID,
                     cMMb = ranLoc_cMMb)
#ranLocsGR <- ranLocsGR[ID_indices]

# Divide features into quantiles based on decreasing nucleotideDiversity
featuresDF <- data.frame(featuresGR,
                         quantile = as.character(""),
                         stringsAsFactors = F)
#featuresDF$nucleotideDiversity[which(is.na(featuresDF$nucleotideDiversity))] <- 0
quantilesStats <- data.frame()
for(k in 1:quantiles) {
  if(k < quantiles) {
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
    featuresDF[ !is.na(featuresDF$nucleotideDiversity) &
                percent_rank(featuresDF$nucleotideDiversity) <= 1-((k-1)/quantiles) &
                percent_rank(featuresDF$nucleotideDiversity) >  1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of features
    featuresDF[ !is.na(featuresDF$nucleotideDiversity) &
                percent_rank(featuresDF$nucleotideDiversity) <= 1-((k-1)/quantiles) &
                percent_rank(featuresDF$nucleotideDiversity) >= 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  }
  write.table(featuresDF[featuresDF$quantile == paste0("Quantile ", k),],
              file = paste0(outDir,
                            "quantile", k, "_of_", quantiles,
                            "_by_nucleotideDiversity_in_",
                            region, "_of_",
                            substring(featureName[1][1], first = 1, last = 5), "_in_",
                            paste0(substring(featureName, first = 10, last = 16),
                                   collapse = "_"), "_",
                            substring(featureName[1][1], first = 18), ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  stats <- data.frame(quantile = as.integer(k),
                      n = as.integer(dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]),
                      mean_width = as.integer(round(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                      total_width = as.integer(sum(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                      mean_nucleotideDiversity = as.numeric(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$nucleotideDiversity, na.rm = T)))
  quantilesStats <- rbind(quantilesStats, stats)
}
write.table(quantilesStats,
            file = paste0(outDir,
                          "summary_", quantiles, "quantiles_by_nucleotideDiversity_in_",
                          region, "_of_",
                          substring(featureName[1][1], first = 1, last = 5), "_in_",
                          paste0(substring(featureName, first = 10, last = 16),
                                 collapse = "_"), "_",
                          substring(featureName[1][1], first = 18), ".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(featuresDF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_nucleotideDiversity_in_",
                          region, "_of_",
                          substring(featureName[1][1], first = 1, last = 5), "_in_",
                          paste0(substring(featureName, first = 10, last = 16),
                                 collapse = "_"), "_",
                          substring(featureName[1][1], first = 18), ".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Divide ranLocs into quantiles based on feature quantile indices
ranLocsDF <- data.frame(ranLocsGR,
                        random = as.character(""),
                        stringsAsFactors = F)
# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$quantile == paste0("Quantile ", k))
})
for(k in 1:quantiles) {
  ranLocsDF[quantileIndices[[k]],]$random <- paste0("Random ", k)
}
write.table(ranLocsDF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_nucleotideDiversity_in_",
                          region, "_of_",
                          substring(featureName[1][1], first = 1, last = 5), "_in_",
                          paste0(substring(featureName, first = 10, last = 16),
                                 collapse = "_"), "_",
                          substring(featureName[1][1], first = 18), "_ranLocs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

## Order features in each quantile by decreasing nucleotideDiversity levels
## to define "row_order" for heatmaps
#combineRowOrders <- function(quantile_bool_list) {
#  do.call("c", lapply(quantile_bool_list, function(x) {
#    quantile_nucleotideDiversity <- rowMeans(nucleotideDiversity[x,], na.rm = T)
#    quantile_nucleotideDiversity[which(is.na(quantile_nucleotideDiversity))] <- 0
#    which(x)[order(quantile_nucleotideDiversity, decreasing = T)]
#  }))
#}
#row_order <- combineRowOrders(quantile_bool_list =
#  lapply(seq_along(1:quantiles), function(k) { 
#    featuresDF$quantile == paste0("Quantile ", k)
#  })
#)
## Confirm row_order is as would be obtained by alternative method
## Note that this alternative 
#stopifnot(identical(row_order,
#                    order(featuresDF$nucleotideDiversity,
#                          decreasing=T)))
#
## Order feature IDs in each quantile by decreasing nucleotideDiversity levels
## for use in GO term enrichment analysis
#listCombineRowOrders <- function(quantile_bool_list) {
#  do.call(list, lapply(quantile_bool_list, function(x) {
#    quantile_nucleotideDiversity <- rowMeans(nucleotideDiversity[x,], na.rm = T)
#    quantile_nucleotideDiversity[which(is.na(quantile_nucleotideDiversity))] <- 0
#    which(x)[order(quantile_nucleotideDiversity, decreasing = T)]
#  }))
#}
#featureIndicesList <- listCombineRowOrders(quantile_bool_list =
#  lapply(seq_along(1:quantiles), function(k) {
#    featuresDF$quantile == paste0("Quantile ", k)
#  })
#)
#stopifnot(identical(row_order,
#                    do.call(c, lapply(featureIndicesList,
#                                      function(x) x))))
## Alternatively, with original ordering:
### Get feature indices for each quantile
##featureIndicesList <- lapply(seq_along(1:quantiles), function(k) {
##  which(featuresDF$quantile == paste0("Quantile ", k))
##})
#
#featureIDsQuantileList <- lapply(seq_along(1:quantiles), function(k) {
#  sub(pattern = "\\.\\d+", replacement = "",
#      x = as.vector(featuresDF[featureIndicesList[[k]],]$featureID))
#})
#sapply(seq_along(featureIDsQuantileList), function(k) {
#  write.table(featureIDsQuantileList[[k]],
#              file = paste0(outDir,
#                            "featureIDs_quantile", k, "_of_", quantiles,
#                            "_by_nucleotideDiversity_in_",
#                            region, "_of_",
#                            substring(featureName[1][1], first = 1, last = 5), "_in_",
#                            paste0(substring(featureName, first = 10, last = 16),
#                                   collapse = "_"), "_",
#                            substring(featureName[1][1], first = 18), ".txt"),
#              quote = F, row.names = F, col.names = F)
#})

## Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
## and sort by decreasing log2mat1RegionRowMeans
#ChIPNames <- c(
#               "ASY1_CS_Rep1_ChIP",
#               "DMC1_Rep1_ChIP",
#               "H2AZ_Rep1_ChIP",
#               "H3K4me3_Rep1_ChIP",
#               "H3K4me1_Rep1_ChIP_SRR8126618",
#               "H3K27ac_Rep1_ChIP_SRR8126621",
#               "H3K27me3_ChIP_SRR6350666",
#               "H3K9me2_Rep1_ChIP",
#               "H3K27me1_Rep1_ChIP"
#              )
#ChIPNamesDir <- c(
#                  "ASY1_CS",
#                  "DMC1",
#                  "H2AZ",
#                  "H3K4me3",
#                  "H3K4me1",
#                  "H3K27ac",
#                  "H3K27me3",
#                  "H3K9me2",
#                  "H3K27me1"
#                 )
#ChIPNamesPlot <- c(
#                   "ASY1",
#                   "DMC1",
#                   "H2A.Z",
#                   "H3K4me3",
#                   "H3K4me1",
#                   "H3K27ac",
#                   "H3K27me3",
#                   "H3K9me2",
#                   "H3K27me1"
#                  )
#ChIPColours <- c(
#                 "purple4",
#                 "green2",
#                 "dodgerblue",
#                 "forestgreen",
#                 "goldenrod1",
#                 "orange",
#                 "navy",
#                 "magenta3",
#                 "firebrick1"
#                )
#otherNames <- c(
#                "MNase_Rep1",
#                "DNaseI_Rep1_SRR8447247",
#                "WT_RNAseq_Rep1_ERR2402974",
#                "WT_RNAseq_Rep2_ERR2402973",
#                "WT_RNAseq_Rep3_ERR2402972"
#               )
#otherNamesDir <- c(
#                   "MNase",
#                   "DNaseI",
#                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
#                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
#                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci"
#                  )
#otherNamesPlot <- c(
#                    "MNase",
#                    "DNaseI",
#                    "RNA-seq Rep1",
#                    "RNA-seq Rep2",
#                    "RNA-seq Rep3"
#                   )
#otherColours <- c(
#                  "darkcyan",
#                  "purple",
#                  "red4",
#                  "red4",
#                  "red4"
#                 )
#sRNANames <- c(
#               "CS+_2_LIB18613_LDI16228"
#              )
#sRNANamesDir <- c(
#                  "sRNAseq_meiocyte_Martin_Moore"
#                 )
#sRNANamesPlot <- c(
#                   "20-nt sRNAs",
#                   "21-nt sRNAs",
#                   "22-nt sRNAs",
#                   "23-nt sRNAs",
#                   "24-nt sRNAs",
#                   "34-nt sRNAs"
#                  )
#sRNAsizes <- c(
#               "20nt",
#               "21nt",
#               "22nt",
#               "23nt",
#               "24nt",
#               "33nt",
#               "34nt"
#              )
#sRNAColours <- c(
#                 "red",
#                 "blue",
#                 "green2",
#                 "darkorange2",
#                 "purple3",
#                 "darkgreen",
#                 "deeppink"
#                )
#DNAmethNames <- c(
#                  "BSseq_Rep8a_SRR6792678"
#                 )
#DNAmethNamesDir <- c(
#                     "BSseq"
#                    )
#DNAmethContexts <- c(
#                     "CpG",
#                     "CHG",
#                     "CHH"
#                    )
#DNAmethNamesPlot <- c(
#                      "mCG",
#                      "mCHG",
#                      "mCHH"
#                     )
#DNAmethColours <- c(
#                    "navy",
#                    "blue",
#                    "deepskyblue1"
#                   )
#
#ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
#  if(ChIPNames[x] %in% c("H3K4me3_ChIP_SRR6350668",
#                         "H3K27me3_ChIP_SRR6350666",
#                         "H3K36me3_ChIP_SRR6350670",
#                         "H3K9ac_ChIP_SRR6350667",
#                         "CENH3_ChIP_SRR1686799")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
#           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else if(ChIPNames[x] %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
#                                "H3K27ac_Rep1_ChIP_SRR8126621")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
#           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else {
#    paste0("/home/ajt200/analysis/wheat/",
#           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  }
#})
#otherDirs <- sapply(seq_along(otherNames), function(x) {
#  if(otherNames[x] %in% c("MNase_Rep1")) {
#    paste0("/home/ajt200/analysis/wheat/",
#           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else if(otherNames[x] %in% c("DNaseI_Rep1_SRR8447247")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
#           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else if(grepl("RNAseq", otherNames[x])) {
#    paste0("/home/ajt200/analysis/wheat/",
#           otherNamesDir[x], "/snakemake_RNAseq_HISAT2/mapped/geneProfiles_subgenomes/matrices/")
#  } else {
#    stop(paste0("otherNames[", x, "] is not compatible with the specified coverage matrix paths"))
#  }
#})
#sRNADirs <- sapply(seq_along(sRNANames), function(x) {
#  if(sRNANames[x] %in% c("CS+_2_LIB18613_LDI16228")) {
#    paste0("/home/ajt200/analysis/wheat/",
#           sRNANamesDir[x], "/snakemake_sRNAseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else {
#    stop(paste0("sRNANames[", x, "] is not compatible with the specified coverage matrix paths"))
#  }
#})
#DNAmethDirs <- sapply(seq_along(DNAmethNames), function(x) {
#  if(DNAmethNames[x] %in% c("BSseq_Rep8a_SRR6792678")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
#           DNAmethNamesDir[x],
#           "/snakemake_BSseq/coverage/geneProfiles_subgenomes/matrices/")
#  } else {
#    stop(paste0("DNAmethNames[", x, "] is not compatible with the specified coverage matrix paths"))
#  }
#})
#
## ChIP
#ChIPmats <- mclapply(seq_along(ChIPNames), function(x) {
#  lapply(seq_along(featureName), function(y) {
#    as.matrix(read.table(paste0(ChIPDirs[x],
#                                ChIPNames[x],
#                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
#                                featureName[y], "_matrix_bin", binName,
#                                "_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(ChIPNames))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#ChIPmats <- mclapply(seq_along(ChIPmats), function(x) {
#  if(length(featureName) == 3) {
#    do.call(rbind, ChIPmats[[x]])
#  } else {
#    ChIPmats[[x]][[1]]
#  }
#}, mc.cores = length(ChIPmats))
#
## Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
## for each matrix depending on library
#log2ChIPmats <- mclapply(seq_along(ChIPmats), function(x) {
#  if(ChIPNames[x] %in% c(
#                         "ASY1_CS_Rep1_ChIP",
#                         "DMC1_Rep1_ChIP",
#                         "H3K4me3_ChIP_SRR6350668",
#                         "H3K27me3_ChIP_SRR6350666",
#                         "H3K36me3_ChIP_SRR6350670",
#                         "H3K9ac_ChIP_SRR6350667",
#                         "H3K4me1_Rep1_ChIP_SRR8126618",
#                         "H3K27ac_Rep1_ChIP_SRR8126621"
#                        )) {
#    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
#    log2((ChIPmats[[x]]+1)/(controlmats[[1]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
#    log2((ChIPmats[[x]]+1)/(controlmats[[2]]+1))
#  }
#}, mc.cores = length(ChIPmats))
#
#for(x in seq_along(log2ChIPmats)) {
#  attr(log2ChIPmats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(log2ChIPmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(log2ChIPmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(log2ChIPmats[[x]], "extend") = c(upstream, downstream)
#  attr(log2ChIPmats[[x]], "smooth") = FALSE
#  attr(log2ChIPmats[[x]], "signal_name") = ChIPNamesPlot[x]
#  attr(log2ChIPmats[[x]], "target_name") = featureName
#  attr(log2ChIPmats[[x]], "target_is_single_point") = FALSE
#  attr(log2ChIPmats[[x]], "background") = 0
#  attr(log2ChIPmats[[x]], "signal_is_categorical") = FALSE
#  class(log2ChIPmats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
#for(x in seq_along(controlmats)) {
#  attr(controlmats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(controlmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(controlmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(controlmats[[x]], "extend") = c(upstream, downstream)
#  attr(controlmats[[x]], "smooth") = FALSE
#  attr(controlmats[[x]], "signal_name") = controlNamesPlot[x]
#  attr(controlmats[[x]], "target_name") = featureName
#  attr(controlmats[[x]], "target_is_single_point") = FALSE
#  attr(controlmats[[x]], "background") = 0
#  attr(controlmats[[x]], "signal_is_categorical") = FALSE
#  class(controlmats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
## other
#othermats <- mclapply(seq_along(otherNames), function(x) {
#  lapply(seq_along(featureName), function(y) {
#    otherFile <- system(paste0("ls ", otherDirs[x],
#                               otherNames[x],
#                               "_MappedOn_wheat_v1.0*", align, "_sort_norm_",
#                               featureName[y], "_matrix_bin", binName,
#                               "_flank", flankName, ".tab"),
#                        intern = T)
#    as.matrix(read.table(otherFile,
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(otherNames))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#othermats <- mclapply(seq_along(othermats), function(x) {
#  if(length(featureName) == 3) {
#    do.call(rbind, othermats[[x]])
#  } else {
#    othermats[[x]][[1]]
#  }
#}, mc.cores = length(othermats))
#
#for(x in seq_along(othermats)) {
#  attr(othermats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(othermats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(othermats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(othermats[[x]], "extend") = c(upstream, downstream)
#  attr(othermats[[x]], "smooth") = FALSE
#  attr(othermats[[x]], "signal_name") = otherNamesPlot[x]
#  attr(othermats[[x]], "target_name") = featureName
#  attr(othermats[[x]], "target_is_single_point") = FALSE
#  attr(othermats[[x]], "background") = 0
#  attr(othermats[[x]], "signal_is_categorical") = FALSE
#  class(othermats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
## sRNA
#sRNAmats <- mclapply(seq_along(sRNAsizes), function(x) {
#  lapply(seq_along(featureName), function(y) {
#    as.matrix(read.table(paste0(sRNADirs,
#                                sRNANames,
#                                "_MappedOn_wheat_v1.0_", align, "_", sRNAsizes[x], "_sort_norm_",
#                                featureName[y], "_matrix_bin", binName,
#                                "_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(sRNAsizes))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#sRNAmats <- mclapply(seq_along(sRNAmats), function(x) {
#  if(length(featureName) == 3) {
#    do.call(rbind, sRNAmats[[x]])
#  } else {
#    sRNAmats[[x]][[1]]
#  }
#}, mc.cores = length(sRNAmats))
#
#for(x in seq_along(sRNAmats)) {
#  attr(sRNAmats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(sRNAmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(sRNAmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(sRNAmats[[x]], "extend") = c(upstream, downstream)
#  attr(sRNAmats[[x]], "smooth") = FALSE
#  attr(sRNAmats[[x]], "signal_name") = sRNANamesPlot[x]
#  attr(sRNAmats[[x]], "target_name") = featureName
#  attr(sRNAmats[[x]], "target_is_single_point") = FALSE
#  attr(sRNAmats[[x]], "background") = 0
#  attr(sRNAmats[[x]], "signal_is_categorical") = FALSE
#  class(sRNAmats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
## DNAmeth
#DNAmethmats <- mclapply(seq_along(DNAmethContexts), function(x) {
#  lapply(seq_along(featureName), function(y) {
#    as.matrix(read.table(paste0(DNAmethDirs,
#                                DNAmethNames,
#                                "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_", DNAmethContexts[x], "_",
#                                featureName[y], "_matrix_bin", binName,
#                                "_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(DNAmethContexts))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#DNAmethmats <- mclapply(seq_along(DNAmethmats), function(x) {
#  if(length(featureName) == 3) {
#    do.call(rbind, DNAmethmats[[x]])
#  } else {
#    DNAmethmats[[x]][[1]]
#  }
#}, mc.cores = length(DNAmethmats))
#
#for(x in seq_along(DNAmethmats)) {
#  attr(DNAmethmats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(DNAmethmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(DNAmethmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(DNAmethmats[[x]], "extend") = c(upstream, downstream)
#  attr(DNAmethmats[[x]], "smooth") = FALSE
#  attr(DNAmethmats[[x]], "signal_name") = DNAmethNamesPlot[x]
#  attr(DNAmethmats[[x]], "target_name") = featureName
#  attr(DNAmethmats[[x]], "target_is_single_point") = FALSE
#  attr(DNAmethmats[[x]], "background") = "NA"
#  attr(DNAmethmats[[x]], "signal_is_categorical") = FALSE
#  class(DNAmethmats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
#
#if(grepl("genes", featureName)) {
#  featureStartLab <- "TSS"
#  featureEndLab <- "TTS"
#} else {
#  featureStartLab <- "Start"
#  featureEndLab <- "End"
#}
#
## Heatmap plotting function
## Note that for plotting heatmaps for individual datasets in separate PDFs,
## must edit this function - print(EnrichedHeatmap(...))
#featureHeatmap <- function(mat,
#                           col_fun,
#                           colour,
#                           datName) {
#  EnrichedHeatmap(mat = mat,
#                  col = col_fun,
#                  column_title = datName,
#                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = colour,
#                                                                                        lwd = 2),
#                                                                              yaxis_side = "left",
#                                                                              yaxis_facing = "left",
#                                                                              yaxis_gp = gpar(fontsize = 10),
#                                                                              pos_line_gp = gpar(col = "black",
#                                                                                                 lty = 2,
#                                                                                                 lwd = 2))),
#                  top_annotation_height = unit(2, "cm"),
#                  width = unit(6, "cm"),
#                  name = datName,
#                  heatmap_legend_param = list(title = datName,
#                                              title_position = "topcenter",
#                                              title_gp = gpar(font = 2, fontsize = 12),
#                                              legend_direction = "horizontal",
#                                              labels_gp = gpar(fontsize = 10)),
#                  axis_name = c(paste0("-", flankNamePlot),
#                                featureStartLab, featureEndLab,
#                                paste0("+", flankNamePlot)),
#                  axis_name_gp = gpar(fontsize = 12),
#                  border = FALSE,
#                  pos_line_gp = gpar(col = "white", lty = 2, lwd = 2),
#                  # If converting into png with pdfTotiffTopng.sh,
#                  # set use_raster to FALSE
#                  #use_raster = FALSE)
#                  use_raster = TRUE, raster_device = "png", raster_quality = 10)
#}
#
## Define heatmap colours
#rich8to6equal <- c("#0000CB", "#0081FF", "#87CEFA", "#FDEE02", "#FFAB00", "#FF3300")
##quantileColours <- c("darkorange1", "green2", "purple3", "deepskyblue")
##quantileColours <- colorRampPalette(c("red", "blue"))(4)
#quantileColours <- c("red", "purple", "blue", "navy")
#
## Create quantile colour block "heatmap"
#quantileBlockhtmp <- Heatmap(featuresDF$quantile,
#                             col = structure(quantileColours,
#                                             names = paste0("Quantile ", 1:quantiles)),
#                             show_row_names = FALSE, show_heatmap_legend = FALSE,
#                             width = unit(3, "mm"), name = "") 
#quantilecMMbheatmap <- Heatmap(featuresDF$cMMb,
#                               cluster_rows = FALSE,
#                               col = colorRamp2(quantile(featuresDF$cMMb,
#                                                         c(0.60, 0.50, 0.40, 0.30),
#                                                         na.rm = T),
#                                                quantileColours), 
#                               na_col = "grey40",
#                               show_row_names = FALSE, show_heatmap_legend = TRUE,
#                               heatmap_legend_param = list(title = "cM/Mb",
#                                                           title_position = "topcenter",
#                                                           title_gp = gpar(font = 2, fontsize = 12),
#                                                           legend_direction = "horizontal",
#                                                           labels_gp = gpar(fontsize = 10)),
#                               width = unit(3, "cm"), name = "")
#quantilecMMbrowAnno <- rowAnnotation(cMMb = anno_points(featuresDF$cMMb,
#                                                        which = "row",
#                                                        size = unit(1, "mm"),
#                                                        gp = gpar(col = "black"),
#                                                        axis_param = list(at = c(0, 1.5), labels = c("0", "1.5")),
#                                                        width = unit(3, "cm")))
## Plot together
#log2ChIPhtmpList <- mclapply(seq_along(ChIPNames), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(log2ChIPmats[[x]],
#                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = log2ChIPmats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = ChIPNamesPlot[x])
#}, mc.cores = length(log2ChIPmats))
#otherhtmpList <- mclapply(seq_along(otherNames), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(othermats[[x]],
#                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = othermats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = otherNamesPlot[x])
#}, mc.cores = length(othermats))
#controlhtmpList <- mclapply(seq_along(controlNames), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(controlmats[[x]],
#                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = controlmats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = controlNamesPlot[x])
#}, mc.cores = length(controlmats))
##sRNAhtmpList <- mclapply(seq_along(sRNANamesPlot), function(x) {
##  ChIP_col_fun <- colorRamp2(quantile(sRNAmats[[x]],
##                                      c(0.998, 0.9982, 0.9984, 0.9986, 0.9988, 0.999),
##                                      na.rm = T),
##                             rich8to6equal)
##  featureHeatmap(mat = sRNAmats[[x]],
##                 col_fun = ChIP_col_fun,
##                 colour = quantileColours,
##                 datName = sRNANamesPlot[x])
##}, mc.cores = length(sRNAmats))
#DNAmethhtmpList <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(DNAmethmats[[x]],
#                                      c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = DNAmethmats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = DNAmethNamesPlot[x])
#}, mc.cores = length(DNAmethmats))
#
#htmpList <- c(quantileBlockhtmp,
#              quantilecMMbheatmap,
#              log2ChIPhtmpList,
#              controlhtmpList[[1]],
#              otherhtmpList,
##              sRNAhtmpList,
#              DNAmethhtmpList)
#htmps <- NULL
#for(x in 1:length(htmpList)) {
#  htmps <- htmps + htmpList[[x]]
#}
#pdf(paste0(plotDir, "log2ChIPcontrol_around_",
#           substring(featureName[1][1], first = 1, last = 5), "_in_",
#           paste0(substring(featureName, first = 10, last = 16),
#                  collapse = "_"), "_",
#           substring(featureName[1][1], first = 18),
#           "_heatmaps_quantiled_by_nucleotideDiversity_in_", region, ".pdf"),
#    width = 3*length(htmpList),
#    height = 10)
#draw(htmps,
#     split = featuresDF$quantile,
#     row_order = row_order,
#     heatmap_legend_side = "bottom",
#     gap = unit(c(1, 1, rep(14, length(htmpList)-2)), "mm")
#    )
#dev.off()
