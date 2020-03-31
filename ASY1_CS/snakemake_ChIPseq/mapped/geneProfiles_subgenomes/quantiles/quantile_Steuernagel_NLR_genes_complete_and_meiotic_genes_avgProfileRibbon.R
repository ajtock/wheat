#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage profiles with 95% CIs around
# gene quantiles; e.g.,
# clusters_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters/cluster1_of_4_by_log2_ASY1_CS_Rep1_ChIP_control_in_promoters_of_genes_in_Agenome_genomewide.txt

# Usage:
# /applications/R/R-3.5.0/bin/Rscript quantile_Steuernagel_NLR_genes_complete_and_meiotic_genes_avgProfileRibbon.R ASY1_CS_Rep1_ChIP ASY1_CS both 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide' 3500 2000 2kb '2 kb' 20 20bp bodies 4 100kb 1 '0.02,0.96'

#libName <- "ASY1_CS_Rep1_ChIP"
#dirName <- "ASY1_CS"
#align <- "both"
#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))
#bodyLength <- 3500
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "bodies"
#quantiles <- 4
#winName <- "100kb"
#minMarkerDist <- 1
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
#                                        split = ",")))

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
align <- args[3]
featureName <- unlist(strsplit(args[4],
                               split = ","))
bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
flankNamePlot <- args[8]
binSize <- as.numeric(args[9])
binName <- args[10]
region <- args[11]
quantiles <- as.numeric(args[12])
winName <- args[13]
minMarkerDist <-  as.numeric(args[14])
legendPos <- as.numeric(unlist(strsplit(args[15],
                                        split = ",")))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(GenomicRanges)

if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  outDir <- paste0("quantiles_by_", libName, "/")
} else {
  outDir <- paste0("quantiles_by_log2_", libName,
                   "_control_in_", region, "/")
}
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Define plot titles
subcat1NamePlot <- paste0("NLRs in ",
                          sub("_\\w+", "", dirName),
                          " quantiles")
subcat2NamePlot <- paste0("Defense in ",
                         sub("_\\w+", "", dirName),
                         " quantiles")
subcat3NamePlot <- paste0("Meiotic in ",
                          sub("_\\w+", "", dirName),
                          " quantiles")
ranFeatNamePlot <- paste0("Random ",
                          substr(featureName[1], start = 1, stop = 4),
                          " quantiles")
ranLocNamePlot <- "Random locus quantiles"

# Define quantile colours
quantileColours <- c("red", "purple", "blue", "navy")

# Define feature start and end labels for plotting
if(grepl("genes", featureName)) {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Genomic definitions
#chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
#chrs <- chrs[-length(chrs)]

chrs <- paste0(rep("chr", 21), rep(1:7, 3),
               c(rep("A", 7), rep("B", 7), rep("D", 7)))

# Load table of features grouped into quantiles
# by decreasing log2(libName/control)
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  featuresDF <- read.table(paste0(outDir, "WesternEurope/features_", quantiles, "quantiles",
                                  "_by_", libName, "_of_",
                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
                                  paste0(substring(featureName, first = 10, last = 16),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 18), "_WesternEurope.txt"),
                           header = T, sep = "\t")
} else {
  featuresDF <- read.table(paste0(outDir, "features_", quantiles, "quantiles",
                                  "_by_log2_", libName, "_control_in_",
                                  region, "_of_",
                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
                                  paste0(substring(featureName, first = 10, last = 16),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 18), ".txt"),
                           header = T, sep = "\t")
}

# Load features to confirm feature (row) ordering in "featuresDF" is the same
# as in "features" (which was used for generating the coverage matrices)
features <- lapply(seq_along(featureName), function(x) {
  read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                    "IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
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
stopifnot(identical(as.character(featuresDF$featureID),
                    as.character(features$V9)))
featureIDs <- sub(pattern = "\\.\\d+", replacement = "",
                  features$V9)

# Load subcat1
subcat1 <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                             "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLR_genes_complete_representative_mRNA.gff3"),
                      header = F)
# Subset features to only those corresponding to subcat1
features_subcat1 <- features[features$V9 %in% subcat1$V9,]
# Get NLR IDs and their row indices in features
subcat1_IDs <- sub(pattern = "\\.\\d+", replacement = "",
                   features_subcat1$V9)
subcat1_ID_indices <- which(featureIDs %in% subcat1_IDs)
# Load subcat1b (to be used for extracting defense response genes not corresponding to NLRs)
subcat1b <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                              "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLR_genes_representative_mRNA.gff3"),
                       header = F)
# Subset features to only those corresponding to subcat1
features_subcat1b <- features[features$V9 %in% subcat1b$V9,]
# Get NLR IDs
subcat1b_IDs <- sub(pattern = "\\.\\d+", replacement = "",
                    features_subcat1b$V9)
# Load table of refseq v1.0 functional annotations to enable removal of NLRs among defense response genes
library(data.table)
IPSanno <- fread(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_FunctionalAnnotation_v1/",
                        "iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned.TAB"),
                 sep = "\t", data.table = F)
# Replace gene model ID decimal suffix (e.g., ".1")
IPSanno$`Gene-ID` <- sub(pattern = "\\.\\d+", replacement = "",
                         x = IPSanno$`Gene-ID`)
# Replace "1G" with "2G" in gene IDs for consistency with v1.1
IPSanno$`Gene-ID` <- sub(pattern = "1G", replacement = "2G",
                         x = IPSanno$`Gene-ID`)
subcat1c_IDs <- sort(unique(c(IPSanno[grepl("NB-ARC", IPSanno$`Interpro-IDs-(Description)`),]$`Gene-ID`,
                              IPSanno[grepl("TIR", IPSanno$`Interpro-IDs-(Description)`),]$`Gene-ID`)))

# Load subcat2
## Load regions associated with local_adaptation (LARs)
#LARs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/He_Akhunov_2019_NatGenet_1000exomes_SNPs/",
#                          "Table_S10_41588_2019_382_MOESM10_ESM.tsv"),
#                   header = T, sep = "\t", stringsAsFactors = F)
#colnames(LARs) <- c("chr", "start", "end", "num_top_SNPs", "env_vars", "num_env_vars")
#LARsGR <- GRanges(seqnames = LARs$chr,
#                  ranges = IRanges(start = LARs$start,
#                                   end = LARs$end),
#                  strand = "*",
#                  env_vars = LARs$env_vars)
#featuresGR <- GRanges(seqnames = featuresDF$seqnames,
#                      ranges = IRanges(start = featuresDF$start,
#                                       end = featuresDF$end),
#                      strand = featuresDF$strand,
#                      featureID = featuresDF$featureID,
#                      quantile = featuresDF$quantile)
#LARsGR_featuresGR_overlap <- findOverlaps(query = LARsGR,
#                                          subject = featuresGR,
#                                          type = "any",
#                                          select = "all",
#                                          ignore.strand = T)
#subcat2 <- featuresGR[unique(subjectHits(LARsGR_featuresGR_overlap))]
# Load functional annotation in order to extract defense response genes
anno <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/RamirezGonzalez_2018_Science_GO_anno/",
                          "RamirezGonzalez_2018_iwgsc_refseqv1.0_OntologiesForGenes_FunctionalAnnotation_HCgenes_in_Agenome_Bgenome_Dgenome_genomewide_GO_IDs_no_chrUn.tsv"),
                   sep = "\t", stringsAsFactors = F)
colnames(anno) <- c("featureID", "GO")
defense_response_indices <- which(grepl(pattern = "GO:0006952", x = anno$GO))
defense_response_indices <- sort(unique(c(defense_response_indices)))
# Retain only defense response genes
# Get defense response gene IDs and their row indices in features
subcat2_IDs <- anno[defense_response_indices,]$featureID
subcat2_IDs <- subcat2_IDs[!(subcat2_IDs %in% unique(c(subcat1b_IDs, subcat1c_IDs)))]
subcat2_ID_indices <- which(featureIDs %in% subcat2_IDs)

# Load subcat3
# Note: these two sets of meiotic genes share 271 common genes that are assigned to a chromosome
meio1 <- read.table(paste0("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/",
                           "Table_S4_meiotic_GO_genes.tsv"),
                    header = T, stringsAsFactors = F)
meio1 <- as.character(meio1$Gene.ID)

meio2 <- read.table(paste0("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/",
                           "Table_S4_meiotic_gene_orthologs.tsv"),
                    header = T, sep = "\t", stringsAsFactors = F)
meio2 <- as.character(meio2$Gene.ID)

meio <- union(meio1, meio2)
subcat3_IDs <- featureIDs[featureIDs %in% meio]
subcat3_ID_indices <- which(featureIDs %in% subcat3_IDs)

nonIDs <- featureIDs[!(featureIDs %in% c(subcat1_IDs, subcat2_IDs, subcat3_IDs))]
# Function to randomly select feature IDs not present in IDs
ran_nonIDs_select <- function(nonIDsChr, n, replaceBool) {
  sample(x = nonIDsChr,
         size = n,
         replace = replaceBool)
}
# Apply ran_nonIDs_select() function on a per-chromosome basis
# and create growing vector of feature IDs called ran_nonIDs
set.seed(9237452)
set.seed(3710475)
subcat1_ran_nonIDs <- NULL
for(i in 1:length(levels(features$V1))) {
  IDsChr <- c(subcat1_IDs)[grepl(paste0("TraesCS",
                                        sub("chr", "", levels(features$V1))[i]),
                                 c(subcat1_IDs))]
  nonIDsChr <- nonIDs[grepl(paste0("TraesCS",
                                   sub("chr", "", levels(features$V1))[i]),
                            nonIDs)]
  ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                     n = length(IDsChr),
                                     replaceBool = FALSE)
  subcat1_ran_nonIDs <- c(subcat1_ran_nonIDs, ran_nonIDsChr)
}
subcat1_ran_nonID_indices <- which(featureIDs %in% subcat1_ran_nonIDs)
set.seed(9237452)
set.seed(3710475)
subcat2_ran_nonIDs <- NULL
for(i in 1:length(levels(features$V1))) {
  IDsChr <- c(subcat2_IDs)[grepl(paste0("TraesCS",
                                        sub("chr", "", levels(features$V1))[i]),
                                 c(subcat2_IDs))]
  nonIDsChr <- nonIDs[grepl(paste0("TraesCS",
                                   sub("chr", "", levels(features$V1))[i]),
                            nonIDs)]
  if(length(nonIDsChr) >= length(IDsChr)) {
    ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                       n = length(IDsChr),
                                       replaceBool = FALSE)
  } else {
    ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                       n = length(IDsChr),
                                       replaceBool = TRUE)
  }
  subcat2_ran_nonIDs <- c(subcat2_ran_nonIDs, ran_nonIDsChr)
}
subcat2_ran_nonID_indices <- which(featureIDs %in% subcat2_ran_nonIDs)
set.seed(9237452)
set.seed(3710475)
subcat3_ran_nonIDs <- NULL
for(i in 1:length(levels(features$V1))) {
  IDsChr <- c(subcat3_IDs)[grepl(paste0("TraesCS",
                                        sub("chr", "", levels(features$V1))[i]),
                                 c(subcat3_IDs))]
  nonIDsChr <- nonIDs[grepl(paste0("TraesCS",
                                   sub("chr", "", levels(features$V1))[i]),
                            nonIDs)]
  ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                     n = length(IDsChr),
                                     replaceBool = FALSE)
  subcat3_ran_nonIDs <- c(subcat3_ran_nonIDs, ran_nonIDsChr)
}
subcat3_ran_nonID_indices <- which(featureIDs %in% subcat3_ran_nonIDs)


# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresDF$quantile == paste0("Quantile ", k))
})

subcat1_quantileIndices <- lapply(seq_along(quantileIndices), function(k) {
  quantileIndices[[k]][which(quantileIndices[[k]] %in% subcat1_ID_indices)]
})
subcat2_quantileIndices <- lapply(seq_along(quantileIndices), function(k) {
  quantileIndices[[k]][which(quantileIndices[[k]] %in% subcat2_ID_indices)]
})
subcat3_quantileIndices <- lapply(seq_along(quantileIndices), function(k) {
  quantileIndices[[k]][which(quantileIndices[[k]] %in% subcat3_ID_indices)]
})

# Apply ran_nonIDs_select() function on a per-chromosome basis
# and create growing vector of feature IDs called ran_nonIDs
subcat1_randomPCIndices <- lapply(seq_along(subcat1_quantileIndices), function(k) {
  set.seed(9237452)
  set.seed(3710475)
  subcat1_ran_nonIDs <- NULL
  for(i in 1:length(levels(features$V1))) {
    IDsChr <- featureIDs[subcat1_quantileIndices[[k]]][grepl(paste0("TraesCS",
                                                                    sub("chr", "", levels(features$V1))[i]),
                                                             featureIDs[subcat1_quantileIndices[[k]]])]
    nonIDsChr <- nonIDs[grepl(paste0("TraesCS",
                                     sub("chr", "", levels(features$V1))[i]),
                              nonIDs)]
    ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                       n = length(IDsChr),
                                       replaceBool = FALSE)
    subcat1_ran_nonIDs <- c(subcat1_ran_nonIDs, ran_nonIDsChr)
  }
  subcat1_ran_nonID_indices <- which(featureIDs %in% subcat1_ran_nonIDs)
  return(subcat1_ran_nonID_indices)
})
subcat2_randomPCIndices <- lapply(seq_along(subcat2_quantileIndices), function(k) {
  set.seed(9237452)
  set.seed(3710475)
  subcat2_ran_nonIDs <- NULL
  for(i in 1:length(levels(features$V1))) {
    IDsChr <- featureIDs[subcat2_quantileIndices[[k]]][grepl(paste0("TraesCS",
                                                                    sub("chr", "", levels(features$V1))[i]),
                                                             featureIDs[subcat2_quantileIndices[[k]]])]
    nonIDsChr <- nonIDs[grepl(paste0("TraesCS",
                                     sub("chr", "", levels(features$V1))[i]),
                              nonIDs)]
    if(length(nonIDsChr) >= length(IDsChr)) {
      ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                         n = length(IDsChr),
                                         replaceBool = FALSE)
    } else {
      ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                         n = length(IDsChr),
                                         replaceBool = TRUE)
    }
    subcat2_ran_nonIDs <- c(subcat2_ran_nonIDs, ran_nonIDsChr)
  }
  subcat2_ran_nonID_indices <- which(featureIDs %in% subcat2_ran_nonIDs)
  return(subcat2_ran_nonID_indices)
})
subcat3_randomPCIndices <- lapply(seq_along(subcat3_quantileIndices), function(k) {
  set.seed(9237452)
  set.seed(3710475)
  subcat3_ran_nonIDs <- NULL
  for(i in 1:length(levels(features$V1))) {
    IDsChr <- featureIDs[subcat3_quantileIndices[[k]]][grepl(paste0("TraesCS",
                                                                    sub("chr", "", levels(features$V1))[i]),
                                                             featureIDs[subcat3_quantileIndices[[k]]])]
    nonIDsChr <- nonIDs[grepl(paste0("TraesCS",
                                     sub("chr", "", levels(features$V1))[i]),
                              nonIDs)]
    ran_nonIDsChr <- ran_nonIDs_select(nonIDsChr = nonIDsChr,
                                       n = length(IDsChr),
                                       replaceBool = FALSE)
    subcat3_ran_nonIDs <- c(subcat3_ran_nonIDs, ran_nonIDsChr)
  }
  subcat3_ran_nonID_indices <- which(featureIDs %in% subcat3_ran_nonIDs)
  return(subcat3_ran_nonID_indices)
})

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
# and sort by decreasing log2mat1RegionRowMeans
ChIPNames <- c(
               "ASY1_CS_Rep1_ChIP",
               "DMC1_Rep1_ChIP",
               "H2AZ_Rep1_ChIP",
               "H3K4me3_Rep1_ChIP",
               "H3K4me1_Rep1_ChIP_SRR8126618",
               "H3K27ac_Rep1_ChIP_SRR8126621",
               "H3K27me3_ChIP_SRR6350666",
               "H3K9me2_Rep1_ChIP",
               "H3K27me1_Rep1_ChIP"
              )
ChIPNamesDir <- c(
                  "ASY1_CS",
                  "DMC1",
                  "H2AZ",
                  "H3K4me3",
                  "H3K4me1",
                  "H3K27ac",
                  "H3K27me3",
                  "H3K9me2",
                  "H3K27me1"
                 )
log2ChIPNamesPlot <- c(
                       "ASY1",
                       "DMC1",
                       "H2A.Z",
                       "H3K4me3",
                       "H3K4me1",
                       "H3K27ac",
                       "H3K27me3",
                       "H3K9me2",
                       "H3K27me1"
                      )
log2ChIPColours <- c(
                     "purple4",
                     "green2",
                     "dodgerblue",
                     "forestgreen",
                     "goldenrod1",
                     "orange",
                     "navy",
                     "magenta3",
                     "firebrick1"
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
controlColours <- c(
                    "grey40",
                    "darkcyan"
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

otherNames <- c(
                "MNase_Rep1",
                "DNaseI_Rep1_SRR8447247",
                "WT_RNAseq_Rep1_ERR2402974",
                "WT_RNAseq_Rep2_ERR2402973",
                "WT_RNAseq_Rep3_ERR2402972"
               )
otherNamesDir <- c(
                   "MNase",
                   "DNaseI",
                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci"
                  )
otherNamesPlot <- c(
                    "MNase",
                    "DNaseI",
                    "RNA-seq Rep1",
                    "RNA-seq Rep2",
                    "RNA-seq Rep3"
                   )
otherColours <- c(
                  "darkcyan",
                  "purple",
                  "red4",
                  "red4",
                  "red4"
                 )
otherDirs <- sapply(seq_along(otherNames), function(x) {
  if(otherNames[x] %in% c("MNase_Rep1")) {
    paste0("/home/ajt200/analysis/wheat/",
           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else if(otherNames[x] %in% c("DNaseI_Rep1_SRR8447247")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
  } else if(grepl("RNAseq", otherNames[x])) {
    paste0("/home/ajt200/analysis/wheat/",
           otherNamesDir[x], "/snakemake_RNAseq_HISAT2/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("otherNames[", x, "] is not compatible with the specified coverage matrix paths"))
  }
})

sRNANames <- c(
               "CS+_2_LIB18613_LDI16228"
              )
sRNANamesDir <- c(
                  "sRNAseq_meiocyte_Martin_Moore"
                 )
sRNANamesPlot <- c(
                   "20-nt sRNAs",
                   "21-nt sRNAs",
                   "22-nt sRNAs",
                   "23-nt sRNAs",
                   "24-nt sRNAs",
                   "33-nt sRNAs",
                   "34-nt sRNAs"
                  )
sRNAsizes <- c(
               "20nt",
               "21nt",
               "22nt",
               "23nt",
               "24nt",
               "33nt",
               "34nt"
              )
sRNAColours <- c(
                 "red",
                 "blue",
                 "green2",
                 "darkorange2",
                 "purple3",
                 "darkgreen",
                 "deeppink"
                )
sRNADirs <- sapply(seq_along(sRNANames), function(x) {
  if(sRNANames[x] %in% c("CS+_2_LIB18613_LDI16228")) {
    paste0("/home/ajt200/analysis/wheat/",
           sRNANamesDir[x], "/snakemake_sRNAseq/mapped/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("sRNANames[", x, "] is not compatible with the specified coverage matrix paths"))
  }
})

DNAmethNames <- c(
                  "BSseq_Rep8a_SRR6792678"
                 )
DNAmethNamesDir <- c(
                     "BSseq"
                    )
DNAmethContexts <- c(
                     "CpG",
                     "CHG",
                     "CHH"
                    )
DNAmethNamesPlot <- c(
                      "mCG",
                      "mCHG",
                      "mCHH"
                     )
DNAmethColours <- c(
                    "navy",
                    "blue",
                    "deepskyblue1"
                   )
DNAmethDirs <- sapply(seq_along(DNAmethNames), function(x) {
  if(DNAmethNames[x] %in% c("BSseq_Rep8a_SRR6792678")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           DNAmethNamesDir[x],
           "/snakemake_BSseq/coverage/geneProfiles_subgenomes/matrices/")
  } else {
    stop(paste0("DNAmethNames[", x, "] is not compatible with the specified coverage matrix paths"))
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

# ranLoc
control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x],
                                controlNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
control_ranLocMats <- mclapply(seq_along(control_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, control_ranLocMats[[x]])
  } else {
    control_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_ranLocMats))

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

# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, ChIP_ranLocMats[[x]])
  } else {
    ChIP_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_ranLocMats))

# Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
# for each matrix depending on library
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
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
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[1]]+1))
  } else {
    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[2]]+1))
  }
}, mc.cores = length(ChIP_ranLocMats))

# Add column names
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
log2ChIP_mats_quantiles <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][subcat1_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_ranLocMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][subcat2_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_ranLocMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][subcat3_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_featureMats[[x]][subcat3_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_ranLocMats[[x]][subcat3_quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(log2ChIP_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats_quantiles), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(log2ChIP_mats_quantiles[[x]][[y]][[k]]),
                 t(log2ChIP_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(log2ChIP_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]])[1])
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  # feature quantiles
  names(summaryDFfeature_list_log2ChIP[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[3]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_log2ChIP[[x]][[4]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[5]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[6]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_log2ChIP[[x]][[7]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[8]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[9]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_log2ChIP into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_log2ChIP  <- mclapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_log2ChIP[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_log2ChIP))
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  # feature quantiles
  summaryDFfeature_log2ChIP[[x]][[1]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_log2ChIP[[x]][[2]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_log2ChIP[[x]][[3]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[3]]))
  # feature quantiles
  summaryDFfeature_log2ChIP[[x]][[4]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[4]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[4]]))
  # feature random groupings
  summaryDFfeature_log2ChIP[[x]][[5]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[5]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[5]]))
  # random loci groupings
  summaryDFfeature_log2ChIP[[x]][[6]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[6]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[6]]))
  # feature quantiles
  summaryDFfeature_log2ChIP[[x]][[7]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[7]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[7]]))
  # feature random groupings
  summaryDFfeature_log2ChIP[[x]][[8]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[8]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[8]]))
  # random loci groupings
  summaryDFfeature_log2ChIP[[x]][[9]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[9]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[9]]))
}

# Define y-axis limits
ymin_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  min(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[4]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[5]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[6]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[7]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[8]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[9]]$CI_lower))
})
ymax_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  max(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[4]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[5]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[6]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[7]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[8]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[9]]$CI_upper))
})

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
### subcat1
## feature
ggObj1_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranFeat
ggObj2_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranLoc
ggObj3_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

### subcat2
## feature
ggObj4_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranFeat
ggObj5_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranLoc
ggObj6_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

### subcat3
## feature
ggObj7_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[7]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[7]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[7]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[7]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat3NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranFeat
ggObj8_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[8]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[8]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[8]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[8]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

## ranLoc
ggObj9_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[9]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[9]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[9]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[9]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = log2ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(log2ChIPNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_log2ChIP,
                                           ggObj2_combined_log2ChIP,
                                           ggObj3_combined_log2ChIP,
                                           ggObj4_combined_log2ChIP,
                                           ggObj5_combined_log2ChIP,
                                           ggObj6_combined_log2ChIP,
                                           ggObj7_combined_log2ChIP,
                                           ggObj8_combined_log2ChIP,
                                           ggObj9_combined_log2ChIP
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(log2ChIPNamesPlot)),
                                                       (length(c(log2ChIPNamesPlot))+1):(length(c(log2ChIPNamesPlot))*2),
                                                       ((length(c(log2ChIPNamesPlot))*2)+1):(length(c(log2ChIPNamesPlot))*3),
                                                       ((length(c(log2ChIPNamesPlot))*3)+1):(length(c(log2ChIPNamesPlot))*4),
                                                       ((length(c(log2ChIPNamesPlot))*4)+1):(length(c(log2ChIPNamesPlot))*5),
                                                       ((length(c(log2ChIPNamesPlot))*5)+1):(length(c(log2ChIPNamesPlot))*6),
                                                       ((length(c(log2ChIPNamesPlot))*6)+1):(length(c(log2ChIPNamesPlot))*7),
                                                       ((length(c(log2ChIPNamesPlot))*7)+1):(length(c(log2ChIPNamesPlot))*8),
                                                       ((length(c(log2ChIPNamesPlot))*8)+1):(length(c(log2ChIPNamesPlot))*9)
                                                      ))
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir,
                "log2ChIPcontrol_avgProfiles_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(log2ChIPNamesPlot)), width = 63, limitsize = FALSE)
} else {
  ggsave(paste0(plotDir,
                "log2ChIPcontrol_avgProfiles_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(log2ChIPNamesPlot)), width = 63, limitsize = FALSE)
}

#### Free up memory by removing no longer required objects
rm(
   ChIP_featureMats, ChIP_ranLocMats,
   control_featureMats, control_ranLocMats,
   log2ChIP_featureMats, log2ChIP_ranLocMats,
   log2ChIP_mats_quantiles,
   wideDFfeature_list_log2ChIP,
   tidyDFfeature_list_log2ChIP,
   summaryDFfeature_list_log2ChIP,
   summaryDFfeature_log2ChIP
  ) 
gc()
#####



## other
# feature
other_featureMats <- mclapply(seq_along(otherNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    otherFile <- system(paste0("ls ", otherDirs[x],
                               otherNames[x],
                               "_MappedOn_wheat_v1.0*", align, "_sort_norm_",
                               featureName[y], "_matrix_bin", binName,
                               "_flank", flankName, ".tab"),
                        intern = T)
    as.matrix(read.table(otherFile,
                         header = F, skip = 3))
  })
}, mc.cores = length(otherNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
other_featureMats <- mclapply(seq_along(other_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, other_featureMats[[x]])
  } else {
    other_featureMats[[x]][[1]]
  }
}, mc.cores = length(other_featureMats))

# ranLoc
other_ranLocMats <- mclapply(seq_along(otherNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    otherFile <- system(paste0("ls ", otherDirs[x],
                               otherNames[x],
                               "_MappedOn_wheat_v1.0*", align, "_sort_norm_",
                               featureName[y], "_ranLoc_matrix_bin", binName,
                               "_flank", flankName, ".tab"),
                        intern = T)
    as.matrix(read.table(otherFile,
                         header = F, skip = 3))
  })
}, mc.cores = length(otherNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
other_ranLocMats <- mclapply(seq_along(other_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, other_ranLocMats[[x]])
  } else {
    other_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(other_ranLocMats))

# Add column names
for(x in seq_along(other_featureMats)) {
  colnames(other_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(other_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
other_mats_quantiles <- mclapply(seq_along(other_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         other_featureMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         other_featureMats[[x]][subcat1_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         other_ranLocMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         other_featureMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         other_featureMats[[x]][subcat2_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         other_ranLocMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         other_featureMats[[x]][subcat3_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         other_featureMats[[x]][subcat3_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         other_ranLocMats[[x]][subcat3_quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(other_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_other <- mclapply(seq_along(other_mats_quantiles), function(x) {
  lapply(seq_along(other_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(other_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(other_mats_quantiles[[x]][[y]][[k]]),
                 t(other_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(other_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_other  <- mclapply(seq_along(wideDFfeature_list_other), function(x) {
  lapply(seq_along(other_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(other_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_other[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_other))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_other)) {
  for(y in seq_along(other_mats_quantiles[[x]])) {
    for(k in seq_along(other_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_other[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_other[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_other[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_other  <- mclapply(seq_along(tidyDFfeature_list_other), function(x) {
  lapply(seq_along(other_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(other_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_other[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_other[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_other[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_other[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_other[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_other[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_other[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_other))

for(x in seq_along(summaryDFfeature_list_other)) {
  for(y in seq_along(other_mats_quantiles[[x]])) {
    for(k in seq_along(other_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_other[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_other[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_other[[x]][[y]][[k]]$window))
      summaryDFfeature_list_other[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_other[[x]][[y]][[k]])[1])
      summaryDFfeature_list_other[[x]][[y]][[k]]$sem <- summaryDFfeature_list_other[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_other[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_other[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_other[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_other[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_other[[x]][[y]][[k]]$sem
      summaryDFfeature_list_other[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_other[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_other[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_other[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_other)) {
  # feature quantiles
  names(summaryDFfeature_list_other[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_other[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_other[[x]][[3]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_other[[x]][[4]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_other[[x]][[5]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_other[[x]][[6]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_other[[x]][[7]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_other[[x]][[8]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_other[[x]][[9]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_other into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_other  <- mclapply(seq_along(summaryDFfeature_list_other), function(x) {
  lapply(seq_along(other_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_other[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_other))
for(x in seq_along(summaryDFfeature_other)) {
  # feature quantiles
  summaryDFfeature_other[[x]][[1]]$quantile <- factor(summaryDFfeature_other[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_other[[x]][[2]]$quantile <- factor(summaryDFfeature_other[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_other[[x]][[3]]$quantile <- factor(summaryDFfeature_other[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[3]]))
  # feature quantiles
  summaryDFfeature_other[[x]][[4]]$quantile <- factor(summaryDFfeature_other[[x]][[4]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[4]]))
  # feature random groupings
  summaryDFfeature_other[[x]][[5]]$quantile <- factor(summaryDFfeature_other[[x]][[5]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[5]]))
  # random loci groupings
  summaryDFfeature_other[[x]][[6]]$quantile <- factor(summaryDFfeature_other[[x]][[6]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[6]]))
  # feature quantiles
  summaryDFfeature_other[[x]][[7]]$quantile <- factor(summaryDFfeature_other[[x]][[7]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[7]]))
  # feature random groupings
  summaryDFfeature_other[[x]][[8]]$quantile <- factor(summaryDFfeature_other[[x]][[8]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[8]]))
  # random loci groupings
  summaryDFfeature_other[[x]][[9]]$quantile <- factor(summaryDFfeature_other[[x]][[9]]$quantile,
                                                         levels = names(summaryDFfeature_list_other[[x]][[9]]))
}

# Define y-axis limits
ymin_list_other <- lapply(seq_along(summaryDFfeature_other), function(x) {
  min(c(summaryDFfeature_other[[x]][[1]]$CI_lower,
        summaryDFfeature_other[[x]][[2]]$CI_lower,
        summaryDFfeature_other[[x]][[3]]$CI_lower,
        summaryDFfeature_other[[x]][[4]]$CI_lower,
        summaryDFfeature_other[[x]][[5]]$CI_lower,
        summaryDFfeature_other[[x]][[6]]$CI_lower,
        summaryDFfeature_other[[x]][[7]]$CI_lower,
        summaryDFfeature_other[[x]][[8]]$CI_lower,
        summaryDFfeature_other[[x]][[9]]$CI_lower))
})
ymax_list_other <- lapply(seq_along(summaryDFfeature_other), function(x) {
  max(c(summaryDFfeature_other[[x]][[1]]$CI_upper,
        summaryDFfeature_other[[x]][[2]]$CI_upper,
        summaryDFfeature_other[[x]][[3]]$CI_upper,
        summaryDFfeature_other[[x]][[4]]$CI_upper,
        summaryDFfeature_other[[x]][[5]]$CI_upper,
        summaryDFfeature_other[[x]][[6]]$CI_upper,
        summaryDFfeature_other[[x]][[7]]$CI_upper,
        summaryDFfeature_other[[x]][[8]]$CI_upper,
        summaryDFfeature_other[[x]][[9]]$CI_upper))
})

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
### subcat1
## feature
ggObj1_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranFeat
ggObj2_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranLoc
ggObj3_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

### subcat2
## feature
ggObj4_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranFeat
ggObj5_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranLoc
ggObj6_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

### subcat3
## feature
ggObj7_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[7]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[7]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[7]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[7]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat3NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranFeat
ggObj8_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[8]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[8]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[8]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[8]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

## ranLoc
ggObj9_combined_other <- mclapply(seq_along(otherNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_other[[x]][[9]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_other[[x]], ymax_list_other[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_other[[x]][[9]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_other[[x]][[9]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_other[[x]][[9]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = otherNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = otherColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(otherNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_other,
                                           ggObj2_combined_other,
                                           ggObj3_combined_other,
                                           ggObj4_combined_other,
                                           ggObj5_combined_other,
                                           ggObj6_combined_other,
                                           ggObj7_combined_other,
                                           ggObj8_combined_other,
                                           ggObj9_combined_other
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(otherNamesPlot)),
                                                       (length(c(otherNamesPlot))+1):(length(c(otherNamesPlot))*2),
                                                       ((length(c(otherNamesPlot))*2)+1):(length(c(otherNamesPlot))*3),
                                                       ((length(c(otherNamesPlot))*3)+1):(length(c(otherNamesPlot))*4),
                                                       ((length(c(otherNamesPlot))*4)+1):(length(c(otherNamesPlot))*5),
                                                       ((length(c(otherNamesPlot))*5)+1):(length(c(otherNamesPlot))*6),
                                                       ((length(c(otherNamesPlot))*6)+1):(length(c(otherNamesPlot))*7),
                                                       ((length(c(otherNamesPlot))*7)+1):(length(c(otherNamesPlot))*8),
                                                       ((length(c(otherNamesPlot))*8)+1):(length(c(otherNamesPlot))*9)
                                                      ))
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir,
                "other_avgProfiles_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(otherNamesPlot)), width = 63, limitsize = FALSE)
} else {
  ggsave(paste0(plotDir,
                "other_avgProfiles_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(otherNamesPlot)), width = 63, limitsize = FALSE)
}

#### Free up memory by removing no longer required objects
rm(
   other_featureMats, other_ranLocMats,
   other_mats_quantiles,
   wideDFfeature_list_other,
   tidyDFfeature_list_other,
   summaryDFfeature_list_other,
   summaryDFfeature_other
  ) 
gc()
#####


## sRNA
# feature
sRNA_featureMats <- mclapply(seq_along(sRNAsizes), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(sRNADirs,
                                sRNANames,
                                "_MappedOn_wheat_v1.0_", align, "_", sRNAsizes[x], "_sort_norm_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(sRNAsizes))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
sRNA_featureMats <- mclapply(seq_along(sRNA_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, sRNA_featureMats[[x]])
  } else {
    sRNA_featureMats[[x]][[1]]
  }
}, mc.cores = length(sRNA_featureMats))

# ranLoc
sRNA_ranLocMats <- mclapply(seq_along(sRNAsizes), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(sRNADirs,
                                sRNANames,
                                "_MappedOn_wheat_v1.0_", align, "_", sRNAsizes[x], "_sort_norm_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(sRNAsizes))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
sRNA_ranLocMats <- mclapply(seq_along(sRNA_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, sRNA_ranLocMats[[x]])
  } else {
    sRNA_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(sRNA_ranLocMats))

# Add column names
for(x in seq_along(sRNA_featureMats)) {
  colnames(sRNA_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(sRNA_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
sRNA_mats_quantiles <- mclapply(seq_along(sRNA_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         sRNA_featureMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         sRNA_featureMats[[x]][subcat1_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         sRNA_ranLocMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         sRNA_featureMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         sRNA_featureMats[[x]][subcat2_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         sRNA_ranLocMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         sRNA_featureMats[[x]][subcat3_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         sRNA_featureMats[[x]][subcat3_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         sRNA_ranLocMats[[x]][subcat3_quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(sRNA_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_sRNA <- mclapply(seq_along(sRNA_mats_quantiles), function(x) {
  lapply(seq_along(sRNA_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(sRNA_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(sRNA_mats_quantiles[[x]][[y]][[k]]),
                 t(sRNA_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(sRNA_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_sRNA  <- mclapply(seq_along(wideDFfeature_list_sRNA), function(x) {
  lapply(seq_along(sRNA_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(sRNA_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_sRNA[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_sRNA))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_sRNA)) {
  for(y in seq_along(sRNA_mats_quantiles[[x]])) {
    for(k in seq_along(sRNA_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_sRNA[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_sRNA[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_sRNA[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_sRNA  <- mclapply(seq_along(tidyDFfeature_list_sRNA), function(x) {
  lapply(seq_along(sRNA_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(sRNA_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_sRNA[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_sRNA[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_sRNA[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_sRNA[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_sRNA[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_sRNA[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_sRNA[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_sRNA))

for(x in seq_along(summaryDFfeature_list_sRNA)) {
  for(y in seq_along(sRNA_mats_quantiles[[x]])) {
    for(k in seq_along(sRNA_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_sRNA[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_sRNA[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_sRNA[[x]][[y]][[k]]$window))
      summaryDFfeature_list_sRNA[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_sRNA[[x]][[y]][[k]])[1])
      summaryDFfeature_list_sRNA[[x]][[y]][[k]]$sem <- summaryDFfeature_list_sRNA[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_sRNA[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_sRNA[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_sRNA[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_sRNA[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_sRNA[[x]][[y]][[k]]$sem
      summaryDFfeature_list_sRNA[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_sRNA[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_sRNA[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_sRNA[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_sRNA)) {
  # feature quantiles
  names(summaryDFfeature_list_sRNA[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_sRNA[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_sRNA[[x]][[3]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_sRNA[[x]][[4]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_sRNA[[x]][[5]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_sRNA[[x]][[6]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_sRNA[[x]][[7]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_sRNA[[x]][[8]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_sRNA[[x]][[9]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_sRNA into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_sRNA  <- mclapply(seq_along(summaryDFfeature_list_sRNA), function(x) {
  lapply(seq_along(sRNA_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_sRNA[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_sRNA))
for(x in seq_along(summaryDFfeature_sRNA)) {
  # feature quantiles
  summaryDFfeature_sRNA[[x]][[1]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_sRNA[[x]][[2]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_sRNA[[x]][[3]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[3]]))
  # feature quantiles
  summaryDFfeature_sRNA[[x]][[4]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[4]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[4]]))
  # feature random groupings
  summaryDFfeature_sRNA[[x]][[5]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[5]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[5]]))
  # random loci groupings
  summaryDFfeature_sRNA[[x]][[6]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[6]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[6]]))
  # feature quantiles
  summaryDFfeature_sRNA[[x]][[7]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[7]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[7]]))
  # feature random groupings
  summaryDFfeature_sRNA[[x]][[8]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[8]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[8]]))
  # random loci groupings
  summaryDFfeature_sRNA[[x]][[9]]$quantile <- factor(summaryDFfeature_sRNA[[x]][[9]]$quantile,
                                                         levels = names(summaryDFfeature_list_sRNA[[x]][[9]]))
}

# Define y-axis limits
ymin_list_sRNA <- lapply(seq_along(summaryDFfeature_sRNA), function(x) {
  min(c(summaryDFfeature_sRNA[[x]][[1]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[2]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[3]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[4]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[5]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[6]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[7]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[8]]$CI_lower,
        summaryDFfeature_sRNA[[x]][[9]]$CI_lower))
})
ymax_list_sRNA <- lapply(seq_along(summaryDFfeature_sRNA), function(x) {
  max(c(summaryDFfeature_sRNA[[x]][[1]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[2]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[3]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[4]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[5]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[6]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[7]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[8]]$CI_upper,
        summaryDFfeature_sRNA[[x]][[9]]$CI_upper))
})

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
### subcat1
## feature
ggObj1_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranFeat
ggObj2_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranLoc
ggObj3_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

### subcat2
## feature
ggObj4_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranFeat
ggObj5_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranLoc
ggObj6_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

### subcat3
## feature
ggObj7_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[7]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[7]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[7]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[7]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat3NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranFeat
ggObj8_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[8]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[8]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[8]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[8]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

## ranLoc
ggObj9_combined_sRNA <- mclapply(seq_along(sRNANamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_sRNA[[x]][[9]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_sRNA[[x]], ymax_list_sRNA[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_sRNA[[x]][[9]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_sRNA[[x]][[9]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_sRNA[[x]][[9]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = sRNANamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = sRNAColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(sRNANamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_sRNA,
                                           ggObj2_combined_sRNA,
                                           ggObj3_combined_sRNA,
                                           ggObj4_combined_sRNA,
                                           ggObj5_combined_sRNA,
                                           ggObj6_combined_sRNA,
                                           ggObj7_combined_sRNA,
                                           ggObj8_combined_sRNA,
                                           ggObj9_combined_sRNA
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(sRNANamesPlot)),
                                                       (length(c(sRNANamesPlot))+1):(length(c(sRNANamesPlot))*2),
                                                       ((length(c(sRNANamesPlot))*2)+1):(length(c(sRNANamesPlot))*3),
                                                       ((length(c(sRNANamesPlot))*3)+1):(length(c(sRNANamesPlot))*4),
                                                       ((length(c(sRNANamesPlot))*4)+1):(length(c(sRNANamesPlot))*5),
                                                       ((length(c(sRNANamesPlot))*5)+1):(length(c(sRNANamesPlot))*6),
                                                       ((length(c(sRNANamesPlot))*6)+1):(length(c(sRNANamesPlot))*7),
                                                       ((length(c(sRNANamesPlot))*7)+1):(length(c(sRNANamesPlot))*8),
                                                       ((length(c(sRNANamesPlot))*8)+1):(length(c(sRNANamesPlot))*9)
                                                      ))
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir,
                "sRNA_avgProfiles_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(sRNANamesPlot)), width = 63, limitsize = FALSE)
} else {
  ggsave(paste0(plotDir,
                "sRNA_avgProfiles_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(sRNANamesPlot)), width = 63, limitsize = FALSE)
}

#### Free up memory by removing no longer required objects
rm(
   sRNA_featureMats, sRNA_ranLocMats,
   sRNA_mats_quantiles,
   wideDFfeature_list_sRNA,
   tidyDFfeature_list_sRNA,
   summaryDFfeature_list_sRNA,
   summaryDFfeature_sRNA
  ) 
gc()
#####



## DNAmeth
# feature
DNAmeth_featureMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs,
                                DNAmethNames,
                                "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_", DNAmethContexts[x], "_",
                                featureName[y], "_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
DNAmeth_featureMats <- mclapply(seq_along(DNAmeth_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, DNAmeth_featureMats[[x]])
  } else {
    DNAmeth_featureMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_featureMats))

# ranLoc
DNAmeth_ranLocMats <- mclapply(seq_along(DNAmethContexts), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0(DNAmethDirs,
                                DNAmethNames,
                                "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_", DNAmethContexts[x], "_",
                                featureName[y], "_ranLoc_matrix_bin", binName,
                                "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(DNAmethContexts))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
DNAmeth_ranLocMats <- mclapply(seq_along(DNAmeth_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, DNAmeth_ranLocMats[[x]])
  } else {
    DNAmeth_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(DNAmeth_ranLocMats))

# Add column names
for(x in seq_along(DNAmeth_featureMats)) {
  colnames(DNAmeth_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(DNAmeth_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
DNAmeth_mats_quantiles <- mclapply(seq_along(DNAmeth_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         DNAmeth_featureMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         DNAmeth_featureMats[[x]][subcat1_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         DNAmeth_ranLocMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         DNAmeth_featureMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         DNAmeth_featureMats[[x]][subcat2_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         DNAmeth_ranLocMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         DNAmeth_featureMats[[x]][subcat3_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         DNAmeth_featureMats[[x]][subcat3_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         DNAmeth_ranLocMats[[x]][subcat3_quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(DNAmeth_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_DNAmeth <- mclapply(seq_along(DNAmeth_mats_quantiles), function(x) {
  lapply(seq_along(DNAmeth_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(DNAmeth_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(DNAmeth_mats_quantiles[[x]][[y]][[k]]),
                 t(DNAmeth_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(DNAmeth_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_DNAmeth  <- mclapply(seq_along(wideDFfeature_list_DNAmeth), function(x) {
  lapply(seq_along(DNAmeth_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(DNAmeth_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_DNAmeth[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_DNAmeth))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_DNAmeth)) {
  for(y in seq_along(DNAmeth_mats_quantiles[[x]])) {
    for(k in seq_along(DNAmeth_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_DNAmeth[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_DNAmeth  <- mclapply(seq_along(tidyDFfeature_list_DNAmeth), function(x) {
  lapply(seq_along(DNAmeth_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(DNAmeth_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_DNAmeth[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_DNAmeth[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_DNAmeth))

for(x in seq_along(summaryDFfeature_list_DNAmeth)) {
  for(y in seq_along(DNAmeth_mats_quantiles[[x]])) {
    for(k in seq_along(DNAmeth_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_DNAmeth[[x]][[y]][[k]]$window))
      summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_DNAmeth[[x]][[y]][[k]])[1])
      summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$sem <- summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$sem
      summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_DNAmeth[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_DNAmeth)) {
  # feature quantiles
  names(summaryDFfeature_list_DNAmeth[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_DNAmeth[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_DNAmeth[[x]][[3]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_DNAmeth[[x]][[4]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_DNAmeth[[x]][[5]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_DNAmeth[[x]][[6]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_DNAmeth[[x]][[7]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_DNAmeth[[x]][[8]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_DNAmeth[[x]][[9]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_DNAmeth into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_DNAmeth  <- mclapply(seq_along(summaryDFfeature_list_DNAmeth), function(x) {
  lapply(seq_along(DNAmeth_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_DNAmeth[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_DNAmeth))
for(x in seq_along(summaryDFfeature_DNAmeth)) {
  # feature quantiles
  summaryDFfeature_DNAmeth[[x]][[1]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_DNAmeth[[x]][[2]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_DNAmeth[[x]][[3]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[3]]))
  # feature quantiles
  summaryDFfeature_DNAmeth[[x]][[4]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[4]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[4]]))
  # feature random groupings
  summaryDFfeature_DNAmeth[[x]][[5]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[5]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[5]]))
  # random loci groupings
  summaryDFfeature_DNAmeth[[x]][[6]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[6]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[6]]))
  # feature quantiles
  summaryDFfeature_DNAmeth[[x]][[7]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[7]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[7]]))
  # feature random groupings
  summaryDFfeature_DNAmeth[[x]][[8]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[8]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[8]]))
  # random loci groupings
  summaryDFfeature_DNAmeth[[x]][[9]]$quantile <- factor(summaryDFfeature_DNAmeth[[x]][[9]]$quantile,
                                                         levels = names(summaryDFfeature_list_DNAmeth[[x]][[9]]))
}

# Define y-axis limits
ymin_list_DNAmeth <- lapply(seq_along(summaryDFfeature_DNAmeth), function(x) {
  min(c(summaryDFfeature_DNAmeth[[x]][[1]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[2]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[3]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[4]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[5]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[6]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[7]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[8]]$CI_lower,
        summaryDFfeature_DNAmeth[[x]][[9]]$CI_lower))
})
ymax_list_DNAmeth <- lapply(seq_along(summaryDFfeature_DNAmeth), function(x) {
  max(c(summaryDFfeature_DNAmeth[[x]][[1]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[2]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[3]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[4]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[5]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[6]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[7]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[8]]$CI_upper,
        summaryDFfeature_DNAmeth[[x]][[9]]$CI_upper))
})

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
### subcat1
## feature
ggObj1_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranFeat
ggObj2_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranLoc
ggObj3_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

### subcat2
## feature
ggObj4_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranFeat
ggObj5_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranLoc
ggObj6_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

### subcat3
## feature
ggObj7_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[7]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[7]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[7]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[7]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat3NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranFeat
ggObj8_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[8]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[8]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[8]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[8]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

## ranLoc
ggObj9_combined_DNAmeth <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_DNAmeth[[x]][[9]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_DNAmeth[[x]], ymax_list_DNAmeth[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_DNAmeth[[x]][[9]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_DNAmeth[[x]][[9]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_DNAmeth[[x]][[9]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = DNAmethNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = DNAmethColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(DNAmethNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_DNAmeth,
                                           ggObj2_combined_DNAmeth,
                                           ggObj3_combined_DNAmeth,
                                           ggObj4_combined_DNAmeth,
                                           ggObj5_combined_DNAmeth,
                                           ggObj6_combined_DNAmeth,
                                           ggObj7_combined_DNAmeth,
                                           ggObj8_combined_DNAmeth,
                                           ggObj9_combined_DNAmeth
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(DNAmethNamesPlot)),
                                                       (length(c(DNAmethNamesPlot))+1):(length(c(DNAmethNamesPlot))*2),
                                                       ((length(c(DNAmethNamesPlot))*2)+1):(length(c(DNAmethNamesPlot))*3),
                                                       ((length(c(DNAmethNamesPlot))*3)+1):(length(c(DNAmethNamesPlot))*4),
                                                       ((length(c(DNAmethNamesPlot))*4)+1):(length(c(DNAmethNamesPlot))*5),
                                                       ((length(c(DNAmethNamesPlot))*5)+1):(length(c(DNAmethNamesPlot))*6),
                                                       ((length(c(DNAmethNamesPlot))*6)+1):(length(c(DNAmethNamesPlot))*7),
                                                       ((length(c(DNAmethNamesPlot))*7)+1):(length(c(DNAmethNamesPlot))*8),
                                                       ((length(c(DNAmethNamesPlot))*8)+1):(length(c(DNAmethNamesPlot))*9)
                                                      ))
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir,
                "DNAmeth_avgProfiles_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(DNAmethNamesPlot)), width = 63, limitsize = FALSE)
} else {
  ggsave(paste0(plotDir,
                "DNAmeth_avgProfiles_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(DNAmethNamesPlot)), width = 63, limitsize = FALSE)
}

#### Free up memory by removing no longer required objects
rm(
   DNAmeth_featureMats, DNAmeth_ranLocMats,
   DNAmeth_mats_quantiles,
   wideDFfeature_list_DNAmeth,
   tidyDFfeature_list_DNAmeth,
   summaryDFfeature_list_DNAmeth,
   summaryDFfeature_DNAmeth
  ) 
gc()
#####


# varietal SNPclasses
SNPclassNames <- c(
                   "all",
                   "missense_variant",
                   "synonymous_variant",
                   "HIGH",
                   "MODERATE",
                   "LOW",
                   "MODIFIER",
                   "upstream_gene_variant",
                   "downstream_gene_variant",
                   "intron_variant",
                   "intergenic",
                   "transition",
                   "transversion"
                  )
SNPclassNamesPlot <- c(
                       "All SNPs",
                       "Missense SNPs",
                       "Synonymous SNPs",
                       "High-impact SNPs",
                       "Moderate-impact SNPs",
                       "Low-impact SNPs",
                       "Modifier-impact SNPs",
                       "SNPs upstream of a gene",
                       "SNPs downstream of a gene",
                       "Intronic SNPs",
                       "Intergenic SNPs",
                       "Transitions",
                       "Transversions"
                      )

# feature
SNPclass_featureMats <- mclapply(seq_along(SNPclassNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                SNPclassNames[x],
                                "_SNPs_around_", featureName[y],
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(SNPclassNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
SNPclass_featureMats <- mclapply(seq_along(SNPclass_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, SNPclass_featureMats[[x]])
  } else {
    SNPclass_featureMats[[x]][[1]]
  }
}, mc.cores = length(SNPclass_featureMats))

# ranLoc
SNPclass_ranLocMats <- mclapply(seq_along(SNPclassNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                SNPclassNames[x],
                                "_SNPs_around_", featureName[y],
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(SNPclassNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
SNPclass_ranLocMats <- mclapply(seq_along(SNPclass_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, SNPclass_ranLocMats[[x]])
  } else {
    SNPclass_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(SNPclass_ranLocMats))

# Add column names
for(x in seq_along(SNPclass_featureMats)) {
  colnames(SNPclass_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(SNPclass_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
SNPclass_mats_quantiles <- mclapply(seq_along(SNPclass_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][subcat1_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         SNPclass_ranLocMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][subcat2_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         SNPclass_ranLocMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][subcat3_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         SNPclass_featureMats[[x]][subcat3_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         SNPclass_ranLocMats[[x]][subcat3_quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(SNPclass_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_SNPclass <- mclapply(seq_along(SNPclass_mats_quantiles), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(SNPclass_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(SNPclass_mats_quantiles[[x]][[y]][[k]]),
                 t(SNPclass_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(SNPclass_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_SNPclass  <- mclapply(seq_along(wideDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(SNPclass_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_SNPclass[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_SNPclass))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_SNPclass)) {
  for(y in seq_along(SNPclass_mats_quantiles[[x]])) {
    for(k in seq_along(SNPclass_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_SNPclass[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_SNPclass  <- mclapply(seq_along(tidyDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(SNPclass_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_SNPclass[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_SNPclass))

for(x in seq_along(summaryDFfeature_list_SNPclass)) {
  for(y in seq_along(SNPclass_mats_quantiles[[x]])) {
    for(k in seq_along(SNPclass_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_SNPclass[[x]][[y]][[k]]$window))
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_SNPclass[[x]][[y]][[k]])[1])
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sem <- summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sem
      summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_SNPclass[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_SNPclass)) {
  # feature quantiles
  names(summaryDFfeature_list_SNPclass[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_SNPclass[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_SNPclass[[x]][[3]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_SNPclass[[x]][[4]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_SNPclass[[x]][[5]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_SNPclass[[x]][[6]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_SNPclass[[x]][[7]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_SNPclass[[x]][[8]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_SNPclass[[x]][[9]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_SNPclass into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_SNPclass  <- mclapply(seq_along(summaryDFfeature_list_SNPclass), function(x) {
  lapply(seq_along(SNPclass_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_SNPclass[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_SNPclass))
for(x in seq_along(summaryDFfeature_SNPclass)) {
  # feature quantiles
  summaryDFfeature_SNPclass[[x]][[1]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_SNPclass[[x]][[2]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_SNPclass[[x]][[3]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[3]]))
  # feature quantiles
  summaryDFfeature_SNPclass[[x]][[4]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[4]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[4]]))
  # feature random groupings
  summaryDFfeature_SNPclass[[x]][[5]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[5]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[5]]))
  # random loci groupings
  summaryDFfeature_SNPclass[[x]][[6]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[6]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[6]]))
  # feature quantiles
  summaryDFfeature_SNPclass[[x]][[7]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[7]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[7]]))
  # feature random groupings
  summaryDFfeature_SNPclass[[x]][[8]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[8]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[8]]))
  # random loci groupings
  summaryDFfeature_SNPclass[[x]][[9]]$quantile <- factor(summaryDFfeature_SNPclass[[x]][[9]]$quantile,
                                                         levels = names(summaryDFfeature_list_SNPclass[[x]][[9]]))
}

# Define y-axis limits
ymin_list_SNPclass <- lapply(seq_along(summaryDFfeature_SNPclass), function(x) {
  min(c(summaryDFfeature_SNPclass[[x]][[1]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[2]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[3]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[4]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[5]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[6]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[7]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[8]]$CI_lower,
        summaryDFfeature_SNPclass[[x]][[9]]$CI_lower))
})
ymax_list_SNPclass <- lapply(seq_along(summaryDFfeature_SNPclass), function(x) {
  max(c(summaryDFfeature_SNPclass[[x]][[1]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[2]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[3]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[4]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[5]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[6]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[7]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[8]]$CI_upper,
        summaryDFfeature_SNPclass[[x]][[9]]$CI_upper))
})

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
### subcat1
## feature
ggObj1_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
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
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranFeat
ggObj2_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranLoc
ggObj3_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

### subcat2
## feature
ggObj4_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
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
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranFeat
ggObj5_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranLoc
ggObj6_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

### subcat3
## feature
ggObj7_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[7]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[7]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[7]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[7]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
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
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat3NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranFeat
ggObj8_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[8]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[8]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[8]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[8]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

## ranLoc
ggObj9_combined_SNPclass <- mclapply(seq_along(SNPclassNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_SNPclass[[x]][[9]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_SNPclass[[x]], ymax_list_SNPclass[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_SNPclass[[x]][[9]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_SNPclass[[x]][[9]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_SNPclass[[x]][[9]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = SNPclassNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(SNPclassNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_SNPclass,
                                           ggObj2_combined_SNPclass,
                                           ggObj3_combined_SNPclass,
                                           ggObj4_combined_SNPclass,
                                           ggObj5_combined_SNPclass,
                                           ggObj6_combined_SNPclass,
                                           ggObj7_combined_SNPclass,
                                           ggObj8_combined_SNPclass,
                                           ggObj9_combined_SNPclass
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(SNPclassNamesPlot)),
                                                       (length(c(SNPclassNamesPlot))+1):(length(c(SNPclassNamesPlot))*2),
                                                       ((length(c(SNPclassNamesPlot))*2)+1):(length(c(SNPclassNamesPlot))*3),
                                                       ((length(c(SNPclassNamesPlot))*3)+1):(length(c(SNPclassNamesPlot))*4),
                                                       ((length(c(SNPclassNamesPlot))*4)+1):(length(c(SNPclassNamesPlot))*5),
                                                       ((length(c(SNPclassNamesPlot))*5)+1):(length(c(SNPclassNamesPlot))*6),
                                                       ((length(c(SNPclassNamesPlot))*6)+1):(length(c(SNPclassNamesPlot))*7),
                                                       ((length(c(SNPclassNamesPlot))*7)+1):(length(c(SNPclassNamesPlot))*8),
                                                       ((length(c(SNPclassNamesPlot))*8)+1):(length(c(SNPclassNamesPlot))*9)
                                                      ))
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir,
                "SNPclass_avgProfiles_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(SNPclassNamesPlot)), width = 63, limitsize = FALSE)
} else {
  ggsave(paste0(plotDir,
                "SNPclass_avgProfiles_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(SNPclassNamesPlot)), width = 63, limitsize = FALSE)
}

#### Free up memory by removing no longer required objects
rm(
   SNPclass_featureMats, SNPclass_ranLocMats,
   SNPclass_mats_quantiles,
   wideDFfeature_list_SNPclass,
   tidyDFfeature_list_SNPclass,
   summaryDFfeature_list_SNPclass,
   summaryDFfeature_SNPclass
  ) 
gc()
#####


# TE superfams
superfamCodes <- c("RLG",
                   "RLC",
                   "RLX",
                   "RIX",
                   "SIX",
                   "DTC",
                   "DTM",
                   "DTX",
                   "DTH",
                   "DMI",
                   "DTT",
                   "DXX",
                   "DTA",
                   "DHH",
                   "XXX")
superfamNames <- c("Gypsy_LTR",
                   "Copia_LTR",
                   "Unclassified_LTR",
                   "LINE",
                   "SINE",
                   "CACTA",
                   "Mutator",
                   "Unclassified_with_TIRs",
                   "Harbinger",
                   "MITE",
                   "Mariner",
                   "Unclassified_class_2",
                   "hAT",
                   "Helitrons",
                   "Unclassified_repeats")
superfamNamesPlot <- c("Gypsy LTR",
                       "Copia LTR",
                       "Unclassified LTR",
                       "LINE",
                       "SINE",
                       "CACTA",
                       "Mutator",
                       "Unclassified with TIRs",
                       "Harbinger",
                       "MITE",
                       "Mariner",
                       "Unclassified class 2",
                       "hAT",
                       "Helitrons",
                       "Unclassified repeats")

# feature
superfam_featureMats <- mclapply(seq_along(superfamNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                superfamNames[x], "_", superfamCodes[x],
                                "_around_", featureName[y],
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(superfamNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
superfam_featureMats <- mclapply(seq_along(superfam_featureMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, superfam_featureMats[[x]])
  } else {
    superfam_featureMats[[x]][[1]]
  }
}, mc.cores = length(superfam_featureMats))

# ranLoc
superfam_ranLocMats <- mclapply(seq_along(superfamNames), function(x) {
  lapply(seq_along(featureName), function(y) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/",
                                superfamNames[x], "_", superfamCodes[x],
                                "_around_", featureName[y],
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = T))
  })
}, mc.cores = length(superfamNames))
# If features from all 3 subgenomes are to be analysed,
# concatenate the 3 corresponding feature coverage matrices
superfam_ranLocMats <- mclapply(seq_along(superfam_ranLocMats), function(x) {
  if(length(featureName) == 3) {
    do.call(rbind, superfam_ranLocMats[[x]])
  } else {
    superfam_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(superfam_ranLocMats))

# Add column names
for(x in seq_along(superfam_featureMats)) {
  colnames(superfam_featureMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(superfam_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
superfam_mats_quantiles <- mclapply(seq_along(superfam_featureMats), function(x) {
  list(
       # feature quantiles
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][subcat1_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         superfam_ranLocMats[[x]][subcat1_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][subcat2_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         superfam_ranLocMats[[x]][subcat2_quantileIndices[[k]],]
       }),
       # feature quantiles
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][subcat3_quantileIndices[[k]],]
       }),
       # feature random groupings
       lapply(1:quantiles, function(k) {
         superfam_featureMats[[x]][subcat3_randomPCIndices[[k]],]
       }),
       # random loci groupings
       lapply(1:quantiles, function(k) {
         superfam_ranLocMats[[x]][subcat3_quantileIndices[[k]],]
       })
      ) 
}, mc.cores = length(superfam_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_superfam <- mclapply(seq_along(superfam_mats_quantiles), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(superfam_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(superfam_mats_quantiles[[x]][[y]][[k]]),
                 t(superfam_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(superfam_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_superfam  <- mclapply(seq_along(wideDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(superfam_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_superfam[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_superfam))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_superfam)) {
  for(y in seq_along(superfam_mats_quantiles[[x]])) {
    for(k in seq_along(superfam_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_superfam[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_superfam[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_superfam  <- mclapply(seq_along(tidyDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(superfam_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_superfam[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_superfam[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_superfam))

for(x in seq_along(summaryDFfeature_list_superfam)) {
  for(y in seq_along(superfam_mats_quantiles[[x]])) {
    for(k in seq_along(superfam_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_superfam[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_superfam[[x]][[y]][[k]]$window))
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_superfam[[x]][[y]][[k]])[1])
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$sem <- summaryDFfeature_list_superfam[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_superfam[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_superfam[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_superfam[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_superfam[[x]][[y]][[k]]$sem
      summaryDFfeature_list_superfam[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_superfam[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_superfam[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_superfam[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_superfam)) {
  # feature quantiles
  names(summaryDFfeature_list_superfam[[x]][[1]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_superfam[[x]][[2]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_superfam[[x]][[3]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_superfam[[x]][[4]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_superfam[[x]][[5]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_superfam[[x]][[6]]) <- randomPCNames
  # feature quantiles
  names(summaryDFfeature_list_superfam[[x]][[7]]) <- quantileNames
  # feature random groupings
  names(summaryDFfeature_list_superfam[[x]][[8]]) <- randomPCNames
  # random loci groupings
  names(summaryDFfeature_list_superfam[[x]][[9]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_superfam into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_superfam  <- mclapply(seq_along(summaryDFfeature_list_superfam), function(x) {
  lapply(seq_along(superfam_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_superfam[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_superfam))
for(x in seq_along(summaryDFfeature_superfam)) {
  # feature quantiles
  summaryDFfeature_superfam[[x]][[1]]$quantile <- factor(summaryDFfeature_superfam[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[1]]))
  # feature random groupings
  summaryDFfeature_superfam[[x]][[2]]$quantile <- factor(summaryDFfeature_superfam[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[2]]))
  # random loci groupings
  summaryDFfeature_superfam[[x]][[3]]$quantile <- factor(summaryDFfeature_superfam[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[3]]))
  # feature quantiles
  summaryDFfeature_superfam[[x]][[4]]$quantile <- factor(summaryDFfeature_superfam[[x]][[4]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[4]]))
  # feature random groupings
  summaryDFfeature_superfam[[x]][[5]]$quantile <- factor(summaryDFfeature_superfam[[x]][[5]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[5]]))
  # random loci groupings
  summaryDFfeature_superfam[[x]][[6]]$quantile <- factor(summaryDFfeature_superfam[[x]][[6]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[6]]))
  # feature quantiles
  summaryDFfeature_superfam[[x]][[7]]$quantile <- factor(summaryDFfeature_superfam[[x]][[7]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[7]]))
  # feature random groupings
  summaryDFfeature_superfam[[x]][[8]]$quantile <- factor(summaryDFfeature_superfam[[x]][[8]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[8]]))
  # random loci groupings
  summaryDFfeature_superfam[[x]][[9]]$quantile <- factor(summaryDFfeature_superfam[[x]][[9]]$quantile,
                                                         levels = names(summaryDFfeature_list_superfam[[x]][[9]]))
}

# Define y-axis limits
ymin_list_superfam <- lapply(seq_along(summaryDFfeature_superfam), function(x) {
  min(c(summaryDFfeature_superfam[[x]][[1]]$CI_lower,
        summaryDFfeature_superfam[[x]][[2]]$CI_lower,
        summaryDFfeature_superfam[[x]][[3]]$CI_lower,
        summaryDFfeature_superfam[[x]][[4]]$CI_lower,
        summaryDFfeature_superfam[[x]][[5]]$CI_lower,
        summaryDFfeature_superfam[[x]][[6]]$CI_lower,
        summaryDFfeature_superfam[[x]][[7]]$CI_lower,
        summaryDFfeature_superfam[[x]][[8]]$CI_lower,
        summaryDFfeature_superfam[[x]][[9]]$CI_lower))
})
ymax_list_superfam <- lapply(seq_along(summaryDFfeature_superfam), function(x) {
  max(c(summaryDFfeature_superfam[[x]][[1]]$CI_upper,
        summaryDFfeature_superfam[[x]][[2]]$CI_upper,
        summaryDFfeature_superfam[[x]][[3]]$CI_upper,
        summaryDFfeature_superfam[[x]][[4]]$CI_upper,
        summaryDFfeature_superfam[[x]][[5]]$CI_upper,
        summaryDFfeature_superfam[[x]][[6]]$CI_upper,
        summaryDFfeature_superfam[[x]][[7]]$CI_upper,
        summaryDFfeature_superfam[[x]][[8]]$CI_upper,
        summaryDFfeature_superfam[[x]][[9]]$CI_upper))
})

# Define legend labels
legendLabs_feature <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeat <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLoc <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
### subcat1
## feature
ggObj1_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
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
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranFeat
ggObj2_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranLoc
ggObj3_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

### subcat2
## feature
ggObj4_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
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
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranFeat
ggObj5_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranLoc
ggObj6_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

### subcat3
## feature
ggObj7_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[7]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[7]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[7]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[7]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
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
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(subcat3NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranFeat
ggObj8_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[8]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[8]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[8]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[8]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeat[[1]]) +
  annotation_custom(legendLabs_ranFeat[[2]]) +
  annotation_custom(legendLabs_ranFeat[[3]]) +
  annotation_custom(legendLabs_ranFeat[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranFeatNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

## ranLoc
ggObj9_combined_superfam <- mclapply(seq_along(superfamNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_superfam[[x]][[9]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_superfam[[x]], ymax_list_superfam[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_superfam[[x]][[9]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_superfam[[x]][[9]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_superfam[[x]][[9]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = superfamNamesPlot[x]) +
  annotation_custom(legendLabs_ranLoc[[1]]) +
  annotation_custom(legendLabs_ranLoc[[2]]) +
  annotation_custom(legendLabs_ranLoc[[3]]) +
  annotation_custom(legendLabs_ranLoc[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 1.0, size = 30)) +
  ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(superfamNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_superfam,
                                           ggObj2_combined_superfam,
                                           ggObj3_combined_superfam,
                                           ggObj4_combined_superfam,
                                           ggObj5_combined_superfam,
                                           ggObj6_combined_superfam,
                                           ggObj7_combined_superfam,
                                           ggObj8_combined_superfam,
                                           ggObj9_combined_superfam
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(superfamNamesPlot)),
                                                       (length(c(superfamNamesPlot))+1):(length(c(superfamNamesPlot))*2),
                                                       ((length(c(superfamNamesPlot))*2)+1):(length(c(superfamNamesPlot))*3),
                                                       ((length(c(superfamNamesPlot))*3)+1):(length(c(superfamNamesPlot))*4),
                                                       ((length(c(superfamNamesPlot))*4)+1):(length(c(superfamNamesPlot))*5),
                                                       ((length(c(superfamNamesPlot))*5)+1):(length(c(superfamNamesPlot))*6),
                                                       ((length(c(superfamNamesPlot))*6)+1):(length(c(superfamNamesPlot))*7),
                                                       ((length(c(superfamNamesPlot))*7)+1):(length(c(superfamNamesPlot))*8),
                                                       ((length(c(superfamNamesPlot))*8)+1):(length(c(superfamNamesPlot))*9)
                                                      ))
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir,
                "TEsuperfam_avgProfiles_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(superfamNamesPlot)), width = 63, limitsize = FALSE)
} else {
  ggsave(paste0(plotDir,
                "TEsuperfam_avgProfiles_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(superfamNamesPlot)), width = 63, limitsize = FALSE)
}

#### Free up memory by removing no longer required objects
rm(
   superfam_featureMats, superfam_ranLocMats,
   superfam_mats_quantiles,
   wideDFfeature_list_superfam,
   tidyDFfeature_list_superfam,
   summaryDFfeature_list_superfam,
   summaryDFfeature_superfam
  ) 
gc()
#####

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_log2ChIP,
                                           ggObj1_combined_other,
                                           ggObj1_combined_sRNA,
                                           ggObj1_combined_DNAmeth,
                                           ggObj1_combined_SNPclass,
                                           ggObj1_combined_superfam,
                                           ggObj2_combined_log2ChIP,
                                           ggObj2_combined_other,
                                           ggObj2_combined_sRNA,
                                           ggObj2_combined_DNAmeth,
                                           ggObj2_combined_SNPclass,
                                           ggObj2_combined_superfam,
                                           ggObj3_combined_log2ChIP,
                                           ggObj3_combined_other,
                                           ggObj3_combined_sRNA,
                                           ggObj3_combined_DNAmeth,
                                           ggObj3_combined_SNPclass,
                                           ggObj3_combined_superfam,
                                           ggObj4_combined_log2ChIP,
                                           ggObj4_combined_other,
                                           ggObj4_combined_sRNA,
                                           ggObj4_combined_DNAmeth,
                                           ggObj4_combined_SNPclass,
                                           ggObj4_combined_superfam,
                                           ggObj5_combined_log2ChIP,
                                           ggObj5_combined_other,
                                           ggObj5_combined_sRNA,
                                           ggObj5_combined_DNAmeth,
                                           ggObj5_combined_SNPclass,
                                           ggObj5_combined_superfam,
                                           ggObj6_combined_log2ChIP,
                                           ggObj6_combined_other,
                                           ggObj6_combined_sRNA,
                                           ggObj6_combined_DNAmeth,
                                           ggObj6_combined_SNPclass,
                                           ggObj6_combined_superfam,
                                           ggObj7_combined_log2ChIP,
                                           ggObj7_combined_other,
                                           ggObj7_combined_sRNA,
                                           ggObj7_combined_DNAmeth,
                                           ggObj7_combined_SNPclass,
                                           ggObj7_combined_superfam,
                                           ggObj8_combined_log2ChIP,
                                           ggObj8_combined_other,
                                           ggObj8_combined_sRNA,
                                           ggObj8_combined_DNAmeth,
                                           ggObj8_combined_SNPclass,
                                           ggObj8_combined_superfam,
                                           ggObj9_combined_log2ChIP,
                                           ggObj9_combined_other,
                                           ggObj9_combined_sRNA,
                                           ggObj9_combined_DNAmeth,
                                           ggObj9_combined_SNPclass,
                                           ggObj9_combined_superfam
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot)),
                                                       (length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*2),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*2)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*3),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*3)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*4),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*4)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*5),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*5)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*6),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*6)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*7),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*7)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*8),
                                                       ((length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*8)+1):(length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot))*9)
                                                      ))
if(libName %in% c("cMMb", "HudsonRM_all", "HudsonRM_syn", "HudsonRM_nonsyn")) {
  ggsave(paste0(plotDir,
                "combined_avgProfiles_around_", quantiles, "quantiles",
                "_by_", libName, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot)), width = 63, limitsize = FALSE)
} else {
  ggsave(paste0(plotDir,
                "combined_avgProfiles_around_", quantiles, "quantiles",
                "_by_log2_", libName, "_control_in_", region, "_of_",
                substring(featureName[1][1], first = 1, last = 5), "_in_",
                paste0(substring(featureName, first = 10, last = 16),
                       collapse = "_"), "_",
                substring(featureName[1][1], first = 18), "_NLR_Defense_Meiotic_genes_v300320.pdf"),
         plot = ggObjGA_combined,
         height = 6.5*length(c(log2ChIPNamesPlot, otherNamesPlot, sRNANamesPlot, DNAmethNamesPlot, SNPclassNamesPlot, superfamNamesPlot)), width = 63, limitsize = FALSE)
}
