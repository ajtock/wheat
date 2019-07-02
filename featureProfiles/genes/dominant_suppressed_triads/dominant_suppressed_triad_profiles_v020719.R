#!/applications/R/R-3.5.0/bin/Rscript

# Plot heatmaps of features sorted by coverage levels in promoters,
# within feature bodies, or in terminators

# Usage:
# /applications/R/R-3.5.0/bin/Rscript dominant_suppressed_triad_profiles_v020719.R 2000 2kb '2 kb' 20 20bp

#flankSize <- 2000 
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- args[2]
flankNamePlot <- args[3]
binSize <- as.numeric(args[4])
binName <- args[5]

# Source functions
source("/projects/ajt200/Rfunctions/wheatPlotFunctions.R")
library(rtracklayer)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

plotDir <- "./plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load matrix of dominant and suppressed gene triads
triads <- read.csv("/home/ajt200/analysis/wheat/H3K27me3/triad_category_per_gene_from_three_leaf.csv",
                   header = T)
print(levels(triads$dataset))
# [1] "HC_CS_no_stress" "HC_Development"

# For the "HC_CS_no_stress" conditions, extract and combine dominant dominant
# and dominant nondominant gene homoeologues within dominant gene triads
noStress_dominant  <- triads[triads$dataset == "HC_CS_no_stress"
                           & triads$general_description == "Dominant",]
noStress_dominant_Adominant  <- triads[triads$dataset == "HC_CS_no_stress"
                                     & triads$general_description == "Dominant"
                                     & grepl("A.dominant", triads$description)
                                     & triads$chr_group == "A",]
noStress_dominant_Bdominant  <- triads[triads$dataset == "HC_CS_no_stress"
                                     & triads$general_description == "Dominant"
                                     & grepl("B.dominant", triads$description)
                                     & triads$chr_group == "B",]
noStress_dominant_Ddominant  <- triads[triads$dataset == "HC_CS_no_stress"
                                     & triads$general_description == "Dominant"
                                     & grepl("D.dominant", triads$description)
                                     & triads$chr_group == "D",]
noStress_dominant_dominant <- rbind(noStress_dominant_Adominant,
                                    noStress_dominant_Bdominant,
                                    noStress_dominant_Ddominant)
noStress_dominant_nondominant <- noStress_dominant[!(noStress_dominant$gene %in%
                                                     noStress_dominant_dominant$gene),]

stopifnot(dim(noStress_dominant_dominant)[1] + 
          dim(noStress_dominant_nondominant)[1] ==
          dim(noStress_dominant)[1])

# For the "HC_CS_no_stress" conditions, extract and combine suppressed suppressed
# and suppressed nonsuppressed gene homoeologues within suppressed gene triads
noStress_suppressed  <- triads[triads$dataset == "HC_CS_no_stress"
                             & triads$general_description == "Suppressed",]
noStress_suppressed_Asuppressed  <- triads[triads$dataset == "HC_CS_no_stress"
                                         & triads$general_description == "Suppressed"
                                         & grepl("A.suppressed", triads$description)
                                         & triads$chr_group == "A",]
noStress_suppressed_Bsuppressed  <- triads[triads$dataset == "HC_CS_no_stress"
                                         & triads$general_description == "Suppressed"
                                         & grepl("B.suppressed", triads$description)
                                         & triads$chr_group == "B",]
noStress_suppressed_Dsuppressed  <- triads[triads$dataset == "HC_CS_no_stress"
                                         & triads$general_description == "Suppressed"
                                         & grepl("D.suppressed", triads$description)
                                         & triads$chr_group == "D",]
noStress_suppressed_suppressed <- rbind(noStress_suppressed_Asuppressed,
                                        noStress_suppressed_Bsuppressed,
                                        noStress_suppressed_Dsuppressed)
noStress_suppressed_nonsuppressed <- noStress_suppressed[!(noStress_suppressed$gene %in%
                                                           noStress_suppressed_suppressed$gene),]

stopifnot(dim(noStress_suppressed_suppressed)[1] + 
          dim(noStress_suppressed_nonsuppressed)[1] ==
          dim(noStress_suppressed)[1])

# For the "HC_CS_no_stress" conditions, extract and combine gene homoeologues
# balanced ("Central") gene triads
noStress_balanced  <- triads[triads$dataset == "HC_CS_no_stress"
                           & triads$general_description == "Central",]

# Sanity check: genes in dominant, suppressed and balanced triads sum to total
stopifnot(dim(noStress_dominant)[1] +
          dim(noStress_suppressed)[1] +
          dim(noStress_balanced)[1] ==
          dim(triads[triads$dataset == "HC_CS_no_stress",])[1])


# For the "HC_Development" conditions, extract and combine dominant dominant
# and dominant nondominant gene homoeologues within dominant gene triads
development_dominant  <- triads[triads$dataset == "HC_Development"
                              & triads$general_description == "Dominant",]
development_dominant_Adominant  <- triads[triads$dataset == "HC_Development"
                                        & triads$general_description == "Dominant"
                                        & grepl("A.dominant", triads$description)
                                        & triads$chr_group == "A",]
development_dominant_Bdominant  <- triads[triads$dataset == "HC_Development"
                                        & triads$general_description == "Dominant"
                                        & grepl("B.dominant", triads$description)
                                        & triads$chr_group == "B",]
development_dominant_Ddominant  <- triads[triads$dataset == "HC_Development"
                                        & triads$general_description == "Dominant"
                                        & grepl("D.dominant", triads$description)
                                        & triads$chr_group == "D",]
development_dominant_dominant <- rbind(development_dominant_Adominant,
                                       development_dominant_Bdominant,
                                       development_dominant_Ddominant)
development_dominant_nondominant <- development_dominant[!(development_dominant$gene %in%
                                                           development_dominant_dominant$gene),]
#stopifnot(dim(development_dominant_dominant)[1] + 
#          dim(development_dominant_nondominant)[1] ==
#          dim(development_dominant)[1])
#Error: dim(development_dominant_dominant)[1] + dim(development_dominant_nondominant)[1] ==  .... is not TRUE
# 5030 (development_dominant_dominant + development_dominant_nondominant) vs 5046 (development_dominant)

# For the "HC_Development" conditions, extract and combine suppressed suppressed
# and suppressed nonsuppressed gene homoeologues within suppressed gene triads
development_suppressed  <- triads[triads$dataset == "HC_Development"
                                & triads$general_description == "Suppressed",]
development_suppressed_Asuppressed  <- triads[triads$dataset == "HC_Development"
                                            & triads$general_description == "Suppressed"
                                            & grepl("A.suppressed", triads$description)
                                            & triads$chr_group == "A",]
development_suppressed_Bsuppressed  <- triads[triads$dataset == "HC_Development"
                                            & triads$general_description == "Suppressed"
                                            & grepl("B.suppressed", triads$description)
                                            & triads$chr_group == "B",]
development_suppressed_Dsuppressed  <- triads[triads$dataset == "HC_Development"
                                            & triads$general_description == "Suppressed"
                                            & grepl("D.suppressed", triads$description)
                                            & triads$chr_group == "D",]
development_suppressed_suppressed <- rbind(development_suppressed_Asuppressed,
                                           development_suppressed_Bsuppressed,
                                           development_suppressed_Dsuppressed)
development_suppressed_nonsuppressed <- development_suppressed[!(development_suppressed$gene %in%
                                                                 development_suppressed_suppressed$gene),]
#stopifnot(dim(development_suppressed_suppressed)[1] + 
#          dim(development_suppressed_nonsuppressed)[1] ==
#          dim(development_suppressed)[1])
# Error: dim(development_suppressed_suppressed)[1] + dim(development_suppressed_nonsuppressed)[1] ==  .... is not TRUE
# 16415 (development_suppressed_suppressed + development_suppressed_nonsuppressed) vs 16539 (development_suppressed)

# For the "HC_Development" conditions, extract and combine gene homoeologues
# balanced ("Central") gene triads
development_balanced  <- triads[triads$dataset == "HC_Development"
                              & triads$general_description == "Central",]

# Sanity check: genes in dominant, suppressed and balanced triads sum to total
stopifnot(dim(development_dominant)[1] +
          dim(development_suppressed)[1] +
          dim(development_balanced)[1] ==
          dim(triads[triads$dataset == "HC_Development",])[1])


# Load representative genes table
genes <- readGFF("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA.gff3")
geneIDs <- gsub("\\.[1-9]+", "", genes$group)
geneIDs <- gsub("2G", "1G", geneIDs)

# For HC_CS_no_stress, obtain row indexes of representative genes table for each gene grouping
noStress_dominant_dominant_idx <- which(geneIDs %in% as.character(noStress_dominant_dominant$gene))
noStress_dominant_nondominant_idx <- which(geneIDs %in% as.character(noStress_dominant_nondominant$gene))
noStress_suppressed_suppressed_idx <- which(geneIDs %in% as.character(noStress_suppressed_suppressed$gene))
noStress_suppressed_nonsuppressed_idx <- which(geneIDs %in% as.character(noStress_suppressed_nonsuppressed$gene))
noStress_balanced_idx <- which(geneIDs %in% as.character(noStress_balanced$gene))

# For HC_Development, obtain row indexes of representative genes table for each gene grouping
development_dominant_dominant_idx <- which(geneIDs %in% as.character(development_dominant_dominant$gene))
development_dominant_nondominant_idx <- which(geneIDs %in% as.character(development_dominant_nondominant$gene))
development_suppressed_suppressed_idx <- which(geneIDs %in% as.character(development_suppressed_suppressed$gene))
development_suppressed_nonsuppressed_idx <- which(geneIDs %in% as.character(development_suppressed_nonsuppressed$gene))
development_balanced_idx <- which(geneIDs %in% as.character(development_balanced$gene))


## Load feature matrices for each sRNA size class and extract gene groupings
#sRNAsizes <- c("20", "21", "22", "23",
#               "24", "33", "34")
#
#sRNAmats <- mclapply(seq_along(sRNAsizes), function(x) {
#  as.matrix(read.table(paste0("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/mapped/geneProfiles/matrices/",
#                              "CS+_2_LIB18613_LDI16228_MappedOn_wheat_v1.0_both_",
#                              sRNAsizes[x], "nt_sort_norm_",
#                              "genes_matrix_bin", binName,
#                              "_flank", flankName, ".tab"),
#                       header = F, skip = 3))
#}, mc.cores = length(sRNAsizes))
#
## HC_CS_no_stress
#sRNAmats_noStress_dominant_dominant <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][noStress_dominant_dominant_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_noStress_dominant_nondominant <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][noStress_dominant_nondominant_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_noStress_suppressed_suppressed <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][noStress_suppressed_suppressed_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_noStress_suppressed_nonsuppressed <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][noStress_suppressed_nonsuppressed_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_noStress_balanced <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][noStress_balanced_idx,])
#}, mc.cores = length(sRNAmats))
#
## HC_Development
#sRNAmats_development_dominant_dominant <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][development_dominant_dominant_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_development_dominant_nondominant <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][development_dominant_nondominant_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_development_suppressed_suppressed <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][development_suppressed_suppressed_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_development_suppressed_nonsuppressed <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][development_suppressed_nonsuppressed_idx,])
#}, mc.cores = length(sRNAmats))
#sRNAmats_development_balanced <- mclapply(seq_along(sRNAmats), function(x) {
#  colMeans(sRNAmats[[x]][development_balanced_idx,])
#}, mc.cores = length(sRNAmats))


# Load feature matrices for each chromatin dataset, log2-transform,
# and sort by decreasing mat1RegionRowMeans
ChIPNames1 <- c(
                "H3K9me2_Rep1_ChIP",
                "H3K4me3_Rep1_ChIP",
                "H3K27me1_Rep1_ChIP",
                "H2AZ_Rep1_ChIP",
                "ASY1_CS_Rep1_ChIP",
                "CENH3_ChIP_SRR1686799"
               )
ChIPNames2 <- c(
                "H3K4me3_ChIP_SRR6350668",
                "H3K9ac_ChIP_SRR6350667",
                "H3K27me3_ChIP_SRR6350666",
                "H3K36me3_ChIP_SRR6350670"
               )
controlNames <- c(
                  "MNase_Rep1",
                  "H3_input_SRR6350669"
                 )
ChIPNamesPlot <- c(
                  "H3K9me2",
                  "H3K4me3",
                  "H3K27me1",
                  "H2A.Z",
                  "ASY1",
                  "CENH3",
                  "H3K4me3 (IWGSC)",
                  "H3K9ac",
                  "H3K27me3",
                  "H3K36me3",
                  "MNase"
                 )
ChIPNames <- c(
               "H3K9me2_Rep1_ChIP",
               "H3K4me3_Rep1_ChIP",
               "H3K27me1_Rep1_ChIP",
               "H2AZ_Rep1_ChIP",
               "ASY1_Rep1_ChIP",
               "CENH3_ChIP_SRR1686799",
               "H3K4me3_ChIP_SRR6350668",
               "H3K9ac_ChIP_SRR6350667",
               "H3K27me3_ChIP_SRR6350666",
               "H3K36me3_ChIP_SRR6350670",
               "MNase_Rep1"
              )
ChIPDirs1 <- c(
               "/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/H3K27me1/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/H2AZ/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/CENH3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"
              )
ChIPDirs2 <- c(
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K4me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K9ac/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
               "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K36me3/snakemake_ChIPseq/mapped/geneProfiles/matrices/"
              )
controlDirs <- c(
                 "/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
                 "/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/input/snakemake_ChIPseq/mapped/geneProfiles/matrices/"             
                )

ChIPmats1 <- mclapply(seq_along(ChIPNames1), function(x) {
  as.matrix(read.table(paste0(ChIPDirs1[x],
                              ChIPNames1[x],
                              "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                              "genes_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames1)) 

ChIPmats2 <- mclapply(seq_along(ChIPNames2), function(x) {
  as.matrix(read.table(paste0(ChIPDirs2[x],
                              ChIPNames2[x],
                              "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                              "genes_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames2)) 

controlmats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                              "genes_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames)) 

log2ChIPmats1 <- mclapply(seq_along(ChIPmats1), function(x) {
  log2((ChIPmats1[[x]]+1)/(controlmats[[1]]+1))
}, mc.cores = length(ChIPmats1))
 
log2ChIPmats2 <- mclapply(seq_along(ChIPmats2), function(x) {
  log2((ChIPmats2[[x]]+1)/(controlmats[[2]]+1))
}, mc.cores = length(ChIPmats2))

log2ChIPmats <- c(log2ChIPmats1, log2ChIPmats2, list(controlmats[[1]]))

# HC_CS_no_stress
log2ChIPmats_noStress_dominant_dominant <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][noStress_dominant_dominant_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_noStress_dominant_nondominant <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][noStress_dominant_nondominant_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_noStress_suppressed_suppressed <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][noStress_suppressed_suppressed_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_noStress_suppressed_nonsuppressed <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][noStress_suppressed_nonsuppressed_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_noStress_balanced <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][noStress_balanced_idx,])
}, mc.cores = length(log2ChIPmats))

# HC_Development
log2ChIPmats_development_dominant_dominant <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][development_dominant_dominant_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_development_dominant_nondominant <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][development_dominant_nondominant_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_development_suppressed_suppressed <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][development_suppressed_suppressed_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_development_suppressed_nonsuppressed <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][development_suppressed_nonsuppressed_idx,])
}, mc.cores = length(log2ChIPmats))
log2ChIPmats_development_balanced <- mclapply(seq_along(log2ChIPmats), function(x) {
  colMeans(log2ChIPmats[[x]][development_balanced_idx,])
}, mc.cores = length(log2ChIPmats))

myColours <- c("grey70", "darkcyan", "cyan", "goldenrod", "goldenrod4")

pdf(paste0(plotDir, "dominant_suppressed_triad_profiles_HC_CS_no_stress_v020719.pdf"),
    height = 2.5*length(log2ChIPmats), width = 3)
par(mfrow = c(length(log2ChIPmats), 1))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

sapply(seq_along(log2ChIPmats), function(x) {
  plotAvgCovDomSupBal(xplot = 1:length(log2ChIPmats_noStress_dominant_dominant[[x]]),
                      balanced = log2ChIPmats_noStress_balanced[[x]],
                      dominant_dominant = log2ChIPmats_noStress_dominant_dominant[[x]],
                      dominant_nondominant = log2ChIPmats_noStress_dominant_nondominant[[x]],
                      suppressed_suppressed = log2ChIPmats_noStress_suppressed_suppressed[[x]],
                      suppressed_nonsuppressed = log2ChIPmats_noStress_suppressed_nonsuppressed[[x]],
                      Ylab = ChIPNamesPlot[x], colours = myColours,
                      flankSize = flankSize, binSize = binSize,
                      flankLabL = paste0("-", flankNamePlot),
                      flankLabR = paste0("+", flankNamePlot),
                      featureStartLab = "TSS", featureEndLab = "TTS",
                      legendLoc = "right",
                      legendLabs = c("Balanced", "Dominant", "Non-dom", "Suppressed", "Non-supp"))
})
#sapply(seq_along(sRNAmats), function(x) {
#  plotAvgCovDomSupBal(xplot = 1:length(sRNAmats_noStress_dominant_dominant[[x]]),
#                      balanced = sRNAmats_noStress_balanced[[x]],
#                      dominant_dominant = sRNAmats_noStress_dominant_dominant[[x]],
#                      dominant_nondominant = sRNAmats_noStress_dominant_nondominant[[x]],
#                      suppressed_suppressed = sRNAmats_noStress_suppressed_suppressed[[x]],
#                      suppressed_nonsuppressed = sRNAmats_noStress_suppressed_nonsuppressed[[x]],
#                      Ylab = paste0(sRNAsizes[x], "-nt sRNAs"), colours = myColours,
#                      flankSize = flankSize, binSize = binSize,
#                      flankLabL = paste0("-", flankNamePlot),
#                      flankLabR = paste0("+", flankNamePlot),
#                      featureStartLab = "TSS", featureEndLab = "TTS",
#                      legendLoc = "top",
#                      legendLabs = c("Balanced", "Dominant", "Non-dom", "Suppressed", "Non-supp"))
#})
dev.off()

pdf(paste0(plotDir, "dominant_suppressed_triad_profiles_HC_Development_v020719.pdf"),
    height = 2.5*length(log2ChIPmats), width = 3)
par(mfrow = c(length(log2ChIPmats), 1))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

sapply(seq_along(log2ChIPmats), function(x) {
  plotAvgCovDomSupBal(xplot = 1:length(log2ChIPmats_development_dominant_dominant[[x]]),
                      balanced = log2ChIPmats_development_balanced[[x]],
                      dominant_dominant = log2ChIPmats_development_dominant_dominant[[x]],
                      dominant_nondominant = log2ChIPmats_development_dominant_nondominant[[x]],
                      suppressed_suppressed = log2ChIPmats_development_suppressed_suppressed[[x]],
                      suppressed_nonsuppressed = log2ChIPmats_development_suppressed_nonsuppressed[[x]],
                      Ylab = ChIPNamesPlot[x], colours = myColours,
                      flankSize = flankSize, binSize = binSize,
                      flankLabL = paste0("-", flankNamePlot),
                      flankLabR = paste0("+", flankNamePlot),
                      featureStartLab = "TSS", featureEndLab = "TTS",
                      legendLoc = "right",
                      legendLabs = c("Balanced", "Dominant", "Non-dom", "Suppressed", "Non-supp"))
})
#sapply(seq_along(sRNAmats), function(x) {
#  plotAvgCovDomSupBal(xplot = 1:length(sRNAmats_development_dominant_dominant[[x]]),
#                      balanced = sRNAmats_development_balanced[[x]],
#                      dominant_dominant = sRNAmats_development_dominant_dominant[[x]],
#                      dominant_nondominant = sRNAmats_development_dominant_nondominant[[x]],
#                      suppressed_suppressed = sRNAmats_development_suppressed_suppressed[[x]],
#                      suppressed_nonsuppressed = sRNAmats_development_suppressed_nonsuppressed[[x]],
#                      Ylab = paste0(sRNAsizes[x], "-nt sRNAs"), colours = myColours,
#                      flankSize = flankSize, binSize = binSize,
#                      flankLabL = paste0("-", flankNamePlot),
#                      flankLabR = paste0("+", flankNamePlot),
#                      featureStartLab = "TSS", featureEndLab = "TTS",
#                      legendLoc = "top",
#                      legendLabs = c("Balanced", "Dominant", "Non-dom", "Suppressed", "Non-supp"))
#})
dev.off()
