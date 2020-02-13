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

#outDir <- paste0("quantiles_by_nucleotideDiversity_in_",
#                 region, "/")
#plotDir <- paste0(outDir, "plots/")
#system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
#system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/He_Akhunov_2019_NatGenet_1000exomes_SNPs/"

# Load 811 exomes-derived SNP matrix from He et al. (2019) Nat. Genet.
if(!(file.exists(paste0(inDir, "all.GP08_mm75_het3_publication01142019.vcf.gz.tbi")))) {
  tabix_build(paste0(inDir, "all.GP08_mm75_het3_publication01142019.vcf.gz"),
              sc = as.integer(1), bc = as.integer(2), ec = as.integer(2), meta = "#", lineskip = as.integer(0))
} else {
  print("*.vcf.gz.tbi index file exists")
}
vcf_handle <- vcf_open(paste0(inDir, "all.GP08_mm75_het3_publication01142019.vcf.gz"))
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

# Specify populations based on regional groupings of accessions
pop_name <- c("NorthAfrica",
              "SubSaharanAfrica",
              "WesternEurope",
              "EasternEurope",
              "MiddleEast",
              "FormerSU",
              "CentralAsia",
              "SouthAsia",
              "EastAsia",
              "NorthAmerica",
              "CentralAmerica",
              "SouthAmerica",
              "Oceania")

pop_list <- lapply(seq_along(pop_name), function(x) {
  as.character(read.table(paste0(inDir, pop_name[x], "_accessions.txt"))[[1]])
})

genomeClass_list <- lapply(seq_along(genomeClass_list), function(x) {
  set.populations(genomeClass_list[[x]],
                  pop_list, diploid = T)
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
## Below will extract stats for first population only
#neutrality_stats_df_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
#  data.frame(get.neutrality(genomeClassSplit_list[[x]],
#                            theta = T)[[1]],
#             stringsAsFactors = F)
#})

genomeClassSplit_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
  F_ST.stats(genomeClassSplit_list[[x]],
             detail = T, mode = "nucleotide", FAST = F)
})
## Below will extract stats for first population only
#F_ST_stats_df_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
#  data.frame(get.F_ST(genomeClassSplit_list[[x]],
#             mode = "nucleotide"),
#             stringsAsFactors = F)
#})

# diversity.stats function generates error when applied to genomeClassSplit_list
# consisting of multiple populations with "keep.site.info = T":
# "Error in `rownames<-`(`*tmp*`, value = nam) :
#   attempt to set 'rownames' on an object with no dimensions"
# However, most diversity stats (but not Pi from Nei) are
# available after running F_ST.stats
# If Pi from Nei is needed, run diversity.stats with "keep.site.info = F"
# (HOWEVER: Pi shouldn't be used when "include.unknown = T"
#  is specified in Whop_readVCF function call, as above)
genomeClassSplit_list <- mclapply(seq_along(genomeClassSplit_list), function(x) {
  diversity.stats(genomeClassSplit_list[[x]],
                  pi = T, keep.site.info = F)
}, mc.cores = detectCores(), mc.preschedule = F)
## Below will extract stats for first population only
#diversity_stats_df_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
#  data.frame(get.diversity(genomeClassSplit_list[[x]],
#                           between = F)[[1]],
#             stringsAsFactors = F)
#})

# For each population, combine statistics into one dataframe
popgen_stats_pop_list <- mclapply(seq_along(pop_list), function(w) {
#for(w in seq_along(pop_list)) {
  popgen_stats_chr_list <- lapply(seq_along(genomeClassSplit_list), function(x) {
    data.frame(chr = as.character(chrs[x]),
               start = as.integer(sub(pattern = " - \\d+", replacement = "",
                                      x = genomeClassSplit_list[[x]]@region.names)),
               end = as.integer(sub(pattern = "\\d+ - ", replacement = "",
                                    x = genomeClassSplit_list[[x]]@region.names)),
               width = as.integer(genomeClassSplit_list[[x]]@n.sites),
               ## OR
               #width = as.integer( ( as.integer(sub(pattern = "\\d+ - ", replacement = "",
               #                                     x = genomeClassSplit_list[[x]]@region.names)) -
               #                      as.integer(sub(pattern = " - \\d+", replacement = "",
               #                                     x = genomeClassSplit_list[[x]]@region.names)) ) + 1 ),
               strand = as.character(features[features$V1 == chrs[x],]$V7),
               ID = as.character(features[features$V1 == chrs[x],]$V9),
               nuc.diversity.within = as.numeric(genomeClassSplit_list[[x]]@nuc.diversity.within[,w]),
               hap.diversity.within = as.numeric(genomeClassSplit_list[[x]]@hap.diversity.within[,w]),
               Pi_Nei = as.numeric(genomeClassSplit_list[[x]]@Pi[,w]),
               nuc.F_ST.vs.all = as.numeric(genomeClassSplit_list[[x]]@nuc.F_ST.vs.all[,w]),
               nucleotide.F_ST = as.numeric(genomeClassSplit_list[[x]]@nucleotide.F_ST[,1]),
               n.segregating.sites = as.numeric(genomeClassSplit_list[[x]]@n.segregating.sites[,w]),
               Tajima.D = as.numeric(genomeClassSplit_list[[x]]@Tajima.D[,w]),
               Rozas.R_2 = as.numeric(genomeClassSplit_list[[x]]@Rozas.R_2[,w]),
               Fu.Li.F = as.numeric(genomeClassSplit_list[[x]]@Fu.Li.F[,w]),
               Fu.Li.D = as.numeric(genomeClassSplit_list[[x]]@Fu.Li.D[,w]),
               theta_Tajima = as.numeric(genomeClassSplit_list[[x]]@theta_Tajima[,w]),
               theta_Watterson = as.numeric(genomeClassSplit_list[[x]]@theta_Watterson[,w]),
               theta_Fu.Li = as.numeric(genomeClassSplit_list[[x]]@theta_Fu.Li[,w]),
               theta_Achaz.Tajima = as.numeric(genomeClassSplit_list[[x]]@theta_Achaz.Tajima[,w]),
               theta_Achaz.Watterson = as.numeric(genomeClassSplit_list[[x]]@theta_Achaz.Watterson[,w]),
               stringsAsFactors = F)
  })
  popgen_stats_pop <- do.call(rbind, popgen_stats_chr_list)
  # Sanity check to ensure popgen_stats feature chromosome IDs and start and end coordinates
  # are identical to those in features data.frame
  if(!identical(popgen_stats_pop$chr, as.character(features$V1))) { 
    stop("popgen_stats chromosome IDs are not identical to those in features data.frame!")
  }
  if(!identical(popgen_stats_pop$start, features$V4)) {
    stop("popgen_stats start coordinates are not identical to those in features data.frame!")
  }
  if(!identical(popgen_stats_pop$end, features$V5)) {
    stop("popgen_stats end coordinates are not identical to those in features data.frame!")
  }
  return(popgen_stats_pop)
}, mc.cores = length(pop_list), mc.preschedule = F)

## For when working with a single population:  
#popgen_stats_df_list <- lapply(seq_along(neutrality_stats_df_list), function(x) {
#  data.frame(chr = as.character(chrs[x]),
#             start = as.integer(sub(pattern = " - \\d+", replacement = "",
#                                    x = genomeClassSplit_list[[x]]@region.names)),
#             end = as.integer(sub(pattern = "\\d+ - ", replacement = "",
#                                  x = genomeClassSplit_list[[x]]@region.names)),
#             width = as.integer( ( as.integer(sub(pattern = "\\d+ - ", replacement = "",
#                                                  x = genomeClassSplit_list[[x]]@region.names)) -
#                                   as.integer(sub(pattern = " - \\d+", replacement = "",
#                                                  x = genomeClassSplit_list[[x]]@region.names)) ) + 1 ),
#             ## OR
#             #width = as.integer(genomeClassSplit_list[[x]]@n.sites),
#             strand = as.character(features[features$V1 == chrs[x],]$V7),
#             ID = as.character(features[features$V1 == chrs[x],]$V9),
#             nuc.diversity.within = as.numeric(diversity_stats_df_list[[x]]$nuc.diversity.within),
#             hap.diversity.within = as.numeric(diversity_stats_df_list[[x]]$hap.diversity.within),
#             Pi = as.numeric(diversity_stats_df_list[[x]]$Pi),
#             hap.F_ST.vs.all = as.numeric(diversity_stats_df_list[[x]]$hap.F_ST.vs.all),
#             nuc.F_ST.vs.all = as.numeric(diversity_stats_df_list[[x]]$nuc.F_ST.vs.all),
#             neutrality_stats_df_list[[x]],
#             stringsAsFactors = F,
#             row.names = as.character(1:dim(neutrality_stats_df_list[[x]])[1]))
#})
#popgen_stats <- do.call(rbind, popgen_stats_df_list)
## Sanity check to ensure popgen_stats feature chromosome IDs and start and end coordinates
## are identical to those in features data.frame
#if(!identical(popgen_stats$chr, features$V1)) { 
#  stop("popgen_stats chromosome IDs are not identical to those in features data.frame!")
#}
#if(!identical(popgen_stats$start, features$V4)) {
#  stop("popgen_stats start coordinates are not identical to those in features data.frame!")
#}
#if(!identical(popgen_stats$end, features$V5)) {
#  stop("popgen_stats end coordinates are not identical to those in features data.frame!")
#}
    

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


# Calculate mean log2(ASY1 ChIP/control) and log2(DMC1 ChIP/control) in
# gene promoters, bodies and terminators
# Extract and save feature IDs for each quantile for further analyses
# (e.g., GO enrichment and average + 95% CI profile plotting).

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
# and sort by decreasing log2mat1RegionRowMeans
ChIPNames <- c(
               "ASY1_CS_Rep1_ChIP",
               "DMC1_Rep1_ChIP"
              )
ChIPNamesDir <- c(
                  "ASY1_CS",
                  "DMC1"
                 )
log2ChIPNamesPlot <- c(
                       "ASY1",
                       "DMC1"
                      )
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  if(ChIPNames[x] %in% c("H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "CENH3_ChIP_SRR1686799")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else if(ChIPNames[x] %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                                "H3K27ac_Rep1_ChIP_SRR8126621")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    paste0("/home/ajt200/analysis/wheat/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  }
})

controlNames <- c(
                  "H3_input_SRR6350669",
                  "MNase_Rep1"
                 )
controlNamesDir <- c(
                     "input",
                     "MNase"
                    )
controlNamesPlot <- c(
                      "Input",
                      "MNase"
                     )
controlDirs <- sapply(seq_along(controlNames), function(x) {
  if(controlNames[x] == "H3_input_SRR6350669") {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           controlNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else if(controlNames[x] == "MNase_Rep1") {
    paste0("/home/ajt200/analysis/wheat/",
           controlNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    if(!(controlNames %in% c("H3_input_SRR6350669", "MNase_Rep1"))) {
      stop(paste0("controlNames[", x, "] is neither H3_input_SRR6350669 nor MNase_Rep1"))
    }
  }
})

## control
# feature
control_featureMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x],
                                controlNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
control_featureMats <- mclapply(seq_along(control_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, control_featureMats[[x]])
  } else {
    control_featureMats[[x]][[1]]
  }
}, mc.cores = length(control_featureMats))

## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, ChIP_featureMats[[x]])
  } else {
    ChIP_featureMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_featureMats))

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(ChIPNames[x] %in% c(
                         "ASY1_CS_Rep1_ChIP",
                         "DMC1_Rep1_ChIP",
                         "H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3K4me1_Rep1_ChIP_SRR8126618",
                         "H3K27ac_Rep1_ChIP_SRR8126621"
                        )) {
    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[1]]+1))
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[2]]+1))
  }
}, mc.cores = length(ChIP_featureMats))

log2ChIP_featureMats_promoters <- lapply(seq_along(log2ChIP_featureMats), function(x) {
  log2ChIP_featureMats[[x]][,(((upstream-1000)/binSize)+1):(upstream/binSize)]
})
log2ChIP_featureMats_bodies <- lapply(seq_along(log2ChIP_featureMats), function(x) {
  log2ChIP_featureMats[[x]][,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
})
log2ChIP_featureMats_terminators <- lapply(seq_along(log2ChIP_featureMats), function(x) {
  log2ChIP_featureMats[[x]][,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
})
log2ChIP_featureMats_promotersRowMeans <- lapply(seq_along(log2ChIP_featureMats_promoters), function(x) {
  rowMeans(log2ChIP_featureMats_promoters[[x]], na.rm = T)
})
log2ChIP_featureMats_bodiesRowMeans <- lapply(seq_along(log2ChIP_featureMats_bodies), function(x) {
  rowMeans(log2ChIP_featureMats_bodies[[x]], na.rm = T)
})
log2ChIP_featureMats_terminatorsRowMeans <- lapply(seq_along(log2ChIP_featureMats_terminators), function(x) {
  rowMeans(log2ChIP_featureMats_terminators[[x]], na.rm = T)
})

# Combine feature coordinates and their corresponding population genetics statistics,
# mean cM/Mb values, and intron and exon numbers
featuresGR_pop_list <- mclapply(seq_along(popgen_stats_pop_list), function(x) {
  featuresGR_pop <- GRanges(featuresGR,
                            featureID = featuresGR$featureID,
                            # Normalise by feature width (per 1 kb)
                            nucleotideDiversity = popgen_stats_pop_list[[x]]$nuc.diversity.within/(popgen_stats_pop_list[[x]]$width/1e3),
                            # The PopGenome developers advise that, where include.unknown = TRUE in the
                            # Whop_readVCF() function call, "diversity.stats: pi and haplotype diversity should not be used"
                            haplotypeDiversity = popgen_stats_pop_list[[x]]$hap.diversity.within,
                            # Normalise by feature width (per 1 kb)
                            Pi_Nei = popgen_stats_pop_list[[x]]$Pi_Nei/(popgen_stats_pop_list[[x]]$width/1e3),
                            nucleotideF_STvAll = popgen_stats_pop_list[[x]]$nuc.F_ST.vs.all,
                            nucleotideF_ST = popgen_stats_pop_list[[x]]$nucleotide.F_ST,
                            nSegregatingSites = popgen_stats_pop_list[[x]]$n.segregating.sites/(popgen_stats_pop_list[[x]]$width/1e3),
                            TajimaD = popgen_stats_pop_list[[x]]$Tajima.D,
                            RozasR2 = popgen_stats_pop_list[[x]]$Rozas.R_2,
                            FuLiF = popgen_stats_pop_list[[x]]$Fu.Li.F,
                            FuLiD = popgen_stats_pop_list[[x]]$Fu.Li.D,
                            ASY1_in_promoters = log2ChIP_featureMats_promotersRowMeans[[1]],
                            DMC1_in_promoters = log2ChIP_featureMats_promotersRowMeans[[2]],
                            ASY1_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[1]],
                            DMC1_in_bodies = log2ChIP_featureMats_bodiesRowMeans[[2]],
                            ASY1_in_terminators = log2ChIP_featureMats_terminatorsRowMeans[[1]],
                            DMC1_in_terminators = log2ChIP_featureMats_terminatorsRowMeans[[2]],
                            cMMb = feature_cMMb,
                            exons = exonNoPerGene,
                            introns = intronNoPerGene)
  #featuresGR_pop <- featuresGR_pop[ID_indices]
  
  # Sanity check to ensure featuresGR_pop feauture IDs, chromosome IDs and start and end coordinates
  # are identical to those in features data.frame
  if(!identical(gsub("\\.\\d+", "", as.character(featuresGR_pop$featureID)),
                gsub("\\.\\d+", "", popgen_stats_pop_list[[x]]$ID))) {
    stop("featuresGR_pop feature IDs are not identical to those in popgen_stats_pop_list[[x]] data.frame!")
  }
  if(!identical(as.character(seqnames(featuresGR_pop)), popgen_stats_pop_list[[x]]$chr)) { 
    stop("featuresGR_pop chromosome IDs are not identical to those in popgen_stats_pop_list[[x]] data.frame!")
  }
  if(!identical(start(featuresGR_pop), popgen_stats_pop_list[[x]]$start)) {
    stop("featuresGR_pop start coordinates are not identical to those in popgen_stats_pop_list[[x]] data.frame!")
  }
  if(!identical(end(featuresGR_pop), popgen_stats_pop_list[[x]]$end)) {
    stop("featuresGR_pop end coordinates are not identical to those in popgen_stats_pop_list[[x]] data.frame!")
  }
  return(featuresGR_pop)
}, mc.cores = length(popgen_stats_pop_list), mc.preschedule = F)

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

# Define first set of ordering factors (population genetics statistics)
# to be used for grouping genes into 2 quantiles
# (3 quantiles can result in very uneven split of genes, or poor separation of values)
quantiles <- 2
orderingFactor <- colnames(data.frame(featuresGR_pop_list[[1]]))[7:16]
outDir <- paste0("quantiles_by_", orderingFactor, "_in_", region, "/")
outDir_list <- lapply(seq_along(outDir), function(w) {
  sapply(seq_along(pop_name), function(x) {
    paste0(outDir[w], pop_name[x], "/")
  })
})
plotDir_list <- lapply(seq_along(outDir), function(w) {
  sapply(seq_along(pop_name), function(x) {
    paste0(outDir_list[[w]][x], "plots/")
  })
})
sapply(seq_along(outDir), function(w) {
 system(paste0("[ -d ", outDir[w], " ] || mkdir ", outDir[w]))
})
sapply(seq_along(outDir), function(w) {
  mclapply(seq_along(pop_name), function(x) {
    system(paste0("[ -d ", outDir_list[[w]][x], " ] || mkdir ", outDir_list[[w]][x]))
  }, mc.cores = length(pop_name), mc.preschedule = F)
})
sapply(seq_along(outDir), function(w) {
  mclapply(seq_along(pop_name), function(x) {
    system(paste0("[ -d ", plotDir_list[[w]][x], " ] || mkdir ", plotDir_list[[w]][x]))
  }, mc.cores = length(pop_name), mc.preschedule = F)
})

# For each population, divide features into quantiles based on decreasing orderingFactor
for(x in 1:length(featuresGR_pop_list)) {
  print(pop_name[x])
  featuresDF <- data.frame(featuresGR_pop_list[[x]],
                           quantile = as.character(""),
                           stringsAsFactors = F)
  mclapply(seq_along(orderingFactor), function(w) {
    print(orderingFactor[w])
    #featuresDF[,which(colnames(featuresDF) == orderingFactor[w])][which(is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]))] <- 0
    quantilesStats <- data.frame()
    for(k in 1:quantiles) {
      if(k < quantiles) {
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
        featuresDF[ !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) &
                    percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >  1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
      } else {
      # Final quantile should span 0 to, e.g., 0.25 proportions of features
        featuresDF[ !is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) &
                    percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >= 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
      }
      write.table(featuresDF[featuresDF$quantile == paste0("Quantile ", k),],
                  file = paste0(outDir_list[[w]][x],
                                "quantile", k, "_of_", quantiles,
                                "_by_", orderingFactor[w], "_in_",
                                region, "_of_",
                                substring(featureName[1][1], first = 1, last = 5), "_in_",
                                paste0(substring(featureName, first = 10, last = 16),
                                       collapse = "_"), "_",
                                substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                  quote = FALSE, sep = "\t", row.names = FALSE)
      stats <- data.frame(quantile = as.integer(k),
                          n = as.integer(dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]),
                          mean_width = as.integer(round(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                          total_width = as.integer(sum(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                          mean_nucleotideDiversity = as.numeric(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),][,which(colnames(featuresDF) == orderingFactor[w])], na.rm = T)))
      quantilesStats <- rbind(quantilesStats, stats)
    }
    write.table(quantilesStats,
                file = paste0(outDir_list[[w]][x],
                              "summary_", quantiles, "quantiles_by_", orderingFactor[w], "_in_",
                              region, "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(featuresDF,
                file = paste0(outDir_list[[w]][x],
                              "features_", quantiles, "quantiles",
                              "_by_", orderingFactor[w], "_in_",
                              region, "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
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
                file = paste0(outDir_list[[w]][x],
                              "features_", quantiles, "quantiles",
                              "_by_", orderingFactor[w], "_in_",
                              region, "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], "_ranLocs.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }, mc.cores = length(orderingFactor), mc.preschedule = F)
}


# Define second set of ordering factors (log2(ChIP/input) in gene promoters, bodies and terminators)
# to be used for grouping genes into 4 quantiles
quantiles <- 4
orderingFactor <- colnames(data.frame(featuresGR_pop_list[[1]]))[17:22]
outDir <- paste0("quantiles_by_", orderingFactor, "/")
outDir_list <- lapply(seq_along(outDir), function(w) {
  sapply(seq_along(pop_name), function(x) {
    paste0(outDir[w], pop_name[x], "/")
  })
})
plotDir_list <- lapply(seq_along(outDir), function(w) {
  sapply(seq_along(pop_name), function(x) {
    paste0(outDir_list[[w]][x], "plots/")
  })
})
sapply(seq_along(outDir), function(w) {
 system(paste0("[ -d ", outDir[w], " ] || mkdir ", outDir[w]))
})
sapply(seq_along(outDir), function(w) {
  mclapply(seq_along(pop_name), function(x) {
    system(paste0("[ -d ", outDir_list[[w]][x], " ] || mkdir ", outDir_list[[w]][x]))
  }, mc.cores = length(pop_name), mc.preschedule = F)
})
sapply(seq_along(outDir), function(w) {
  mclapply(seq_along(pop_name), function(x) {
    system(paste0("[ -d ", plotDir_list[[w]][x], " ] || mkdir ", plotDir_list[[w]][x]))
  }, mc.cores = length(pop_name), mc.preschedule = F)
})

# For each population, divide features into quantiles based on decreasing orderingFactor
for(x in 1:length(featuresGR_pop_list)) {
  print(pop_name[x])
  featuresDF <- data.frame(featuresGR_pop_list[[x]],
                           quantile = as.character(""),
                           stringsAsFactors = F)
  mclapply(seq_along(orderingFactor), function(w) {
    print(orderingFactor[w])
    featuresDF[,which(colnames(featuresDF) == orderingFactor[w])][which(is.na(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]))] <- 0
    quantilesStats <- data.frame()
    for(k in 1:quantiles) {
      if(k < quantiles) {
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
        featuresDF[ percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >  1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
      } else {
      # Final quantile should span 0 to, e.g., 0.25 proportions of features
        featuresDF[ percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(featuresDF[,which(colnames(featuresDF) == orderingFactor[w])]) >= 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
      }
      write.table(featuresDF[featuresDF$quantile == paste0("Quantile ", k),],
                  file = paste0(outDir_list[[w]][x],
                                "quantile", k, "_of_", quantiles,
                                "_by_", orderingFactor[w],
                                "_of_",
                                substring(featureName[1][1], first = 1, last = 5), "_in_",
                                paste0(substring(featureName, first = 10, last = 16),
                                       collapse = "_"), "_",
                                substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                  quote = FALSE, sep = "\t", row.names = FALSE)
      stats <- data.frame(quantile = as.integer(k),
                          n = as.integer(dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]),
                          mean_width = as.integer(round(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                          total_width = as.integer(sum(featuresDF[featuresDF$quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                          mean_nucleotideDiversity = as.numeric(mean(featuresDF[featuresDF$quantile == paste0("Quantile ", k),][,which(colnames(featuresDF) == orderingFactor[w])], na.rm = T)))
      quantilesStats <- rbind(quantilesStats, stats)
    }
    write.table(quantilesStats,
                file = paste0(outDir_list[[w]][x],
                              "summary_", quantiles, "quantiles_by_", orderingFactor[w],
                              "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(featuresDF,
                file = paste0(outDir_list[[w]][x],
                              "features_", quantiles, "quantiles",
                              "_by_", orderingFactor[w],
                              "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
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
                file = paste0(outDir_list[[w]][x],
                              "features_", quantiles, "quantiles",
                              "_by_", orderingFactor[w],
                              "_of_",
                              substring(featureName[1][1], first = 1, last = 5), "_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], "_ranLocs.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }, mc.cores = length(orderingFactor), mc.preschedule = F)
}
