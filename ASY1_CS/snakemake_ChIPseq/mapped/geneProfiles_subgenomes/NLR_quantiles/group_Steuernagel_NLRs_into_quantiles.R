#!/applications/R/R-3.5.0/bin/Rscript

# Group NLR-encoding genes identified in Steuernagel et al. (2019) Plant Physiol.
# into quantiles based on decreasing population genetics statistics,
# coverage for various data sets (e.g., DMC1 and ASY1), and cM/Mb

# Usage:
# /applications/R/R-3.5.0/bin/Rscript group_Steuernagel_NLRs_into_quantiles.R 'genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide'

#featureName <- unlist(strsplit("genes_in_Agenome_genomewide,genes_in_Bgenome_genomewide,genes_in_Dgenome_genomewide",
#                               split = ","))

args <- commandArgs(trailingOnly = T)
featureName <- unlist(strsplit(args[1],
                               split = ","))

library(parallel)
library(dplyr)

# Define wheat exome-seq population name
pop_name <- c(
              "Africa",
              "MiddleEast",
              "Asia",
              "FormerSU",
              "EasternEurope",
              "WesternEurope",
              "NorthAmerica",
              "CentralAmerica",
              "SouthAmerica",
              "Oceania"
             )

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]

# Load table of NLRs including NLRAnnotator-generated motif annotation
# Largest ("repres") set of motifs for each NLR
NLR_repres_motifs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                                       "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLRAnnotator/",
                                       "NLRAnnotator_genes_repres_mRNA_repres_motifs.bed"),
                                header = F, stringsAsFactors = F)
colnames(NLR_repres_motifs) <- c("chr", "start0based", "end", "featureID", "motifs", "strand")
NLR_repres_motifs <- NLR_repres_motifs[NLR_repres_motifs$chr != "chrUn",]
# Concatenated ("concat") sets of motifs for each NLR
NLR_concat_motifs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                                       "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLRAnnotator/",
                                       "NLRAnnotator_genes_repres_mRNA_concat_motifs.bed"),
                                header = F, stringsAsFactors = F)
colnames(NLR_concat_motifs) <- c("chr", "start0based", "end", "featureID", "motifs", "strand")
NLR_concat_motifs <- NLR_concat_motifs[NLR_concat_motifs$chr != "chrUn",]
if(length(featureName) == 1) {
  NLR_repres_motifs <- NLR_repres_motifs[grepl(substr(featureName, 10, 10), NLR_repres_motifs$chr),]
  NLR_concat_motifs <- NLR_concat_motifs[grepl(substr(featureName, 10, 10), NLR_concat_motifs$chr),]
  chrs <- chrs[grepl(substr(featureName, 10, 10), chrs)]
}

# Load table of features (one for each population)
featuresNLR_pop_list <- mclapply(seq_along(pop_name), function(x) {
  featuresDF <- read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/quantiles/",
                                  "quantiles_by_cMMb/", pop_name[x], "/features_4quantiles_by_cMMb_of_",
                                  substring(featureName[1][1], first = 1, last = 5), "_in_",
                                  paste0(substring(featureName, first = 10, last = 16),
                                         collapse = "_"), "_",
                                  substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                           header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
  
  # Get rows in featuresDF corresponding to NLRs
  featuresNLR <- featuresDF[featuresDF$featureID %in% NLR_repres_motifs$featureID,]
  featuresNLR <- featuresNLR[order(featuresNLR$seqnames,
                                   featuresNLR$start,
                                   featuresNLR$end),]
  stopifnot(identical(NLR_repres_motifs$featureID, featuresNLR$featureID))
  stopifnot(identical(NLR_concat_motifs$featureID, featuresNLR$featureID))
  featuresNLR$clustered <- ""
  featuresNLR$cluster <- as.integer("")
  featuresNLR$cluster_members <- as.integer("")
  # Determine if each NLR-encoding gene is part of a cluster of NLR-encoding genes
  # (i.e., not interrupted by more than eight genes encoding non-NLRs,
  # as defined in Richly et al. (2002) Mol. Biol. Evol. 19: https://doi.org/10.1093/oxfordjournals.molbev.a003984 )
  # Assign cluster number to each NLR-encoding gene within a cluster
  featuresNLR_cl <- NULL
  for(i in seq_along(chrs)) {
    print(chrs[i])
    featuresNLR_chr <- featuresNLR[featuresNLR$seqnames == chrs[i],]
    # First NLR in a chromosome
    if( as.integer(rownames(featuresNLR_chr[2,])) -
        as.integer(rownames(featuresNLR_chr[1,])) - 1 <=
        8 ) {
      featuresNLR_chr[1,]$clustered <- "yes"
      featuresNLR_chr[1,]$cluster <- as.integer(1)
    } else {
      featuresNLR_chr[1,]$clustered <- "no"
    }
    # All other NLRs in a chromosome, except last 
    for(y in 2:(dim(featuresNLR_chr)[1]-1)) {
      if( as.integer(rownames(featuresNLR_chr[y,])) -
          as.integer(rownames(featuresNLR_chr[y-1,])) - 1 <=
          8 ) {
        featuresNLR_chr[y,]$clustered <- "yes"
        featuresNLR_chr[y,]$cluster <- as.integer(featuresNLR_chr[y-1,]$cluster) 
      } else if ( as.integer(rownames(featuresNLR_chr[y+1,])) -
                  as.integer(rownames(featuresNLR_chr[y,])) - 1 <=
                  8 ) {
        featuresNLR_chr[y,]$clustered <- "yes"
        if( length(featuresNLR_chr[!is.na(featuresNLR_chr$cluster),]$cluster) > 0 ) {
          featuresNLR_chr[y,]$cluster <- as.integer(featuresNLR_chr[!is.na(featuresNLR_chr$cluster),]$cluster[length(featuresNLR_chr[!is.na(featuresNLR_chr$cluster),]$cluster)] + 1)
        } else {
          featuresNLR_chr[y,]$cluster <- as.integer(1)
        }
      } else {
        featuresNLR_chr[y,]$clustered <- "no"
      }
    }
    # Last NLR in a chromosome
    if( as.integer(rownames(featuresNLR_chr[dim(featuresNLR_chr)[1],])) -
        as.integer(rownames(featuresNLR_chr[dim(featuresNLR_chr)[1]-1,])) - 1 <=
        8 ) {
      featuresNLR_chr[dim(featuresNLR_chr)[1],]$clustered <- "yes"
      featuresNLR_chr[dim(featuresNLR_chr)[1],]$cluster <- as.integer(featuresNLR_chr[dim(featuresNLR_chr)[1]-1,]$cluster)
    } else {
      featuresNLR_chr[dim(featuresNLR_chr)[1],]$clustered <- "no"
    }
    for(z in 1:dim(featuresNLR_chr)[1]) {
      if( featuresNLR_chr[z,]$clustered == "yes" ) {
        featuresNLR_chr[z,]$cluster_members <- as.integer(dim(featuresNLR_chr[which(featuresNLR_chr$cluster == featuresNLR_chr[z,]$cluster),])[1])
      }
    }
    featuresNLR_cl <- rbind(featuresNLR_cl, featuresNLR_chr)
  }
  featuresNLR <- featuresNLR_cl
  stopifnot(identical(NLR_repres_motifs$featureID, featuresNLR$featureID))
  stopifnot(identical(NLR_concat_motifs$featureID, featuresNLR$featureID))
  featuresNLR$repres_motifs <- NLR_repres_motifs$motifs
  featuresNLR$concat_motifs <- NLR_concat_motifs$motifs
  return(featuresNLR)
}, mc.cores = length(pop_name), mc.preschedule = F)

# Define set of ordering factors to be used for grouping genes into 4 quantiles
orderingFactor <- colnames(featuresNLR_pop_list[[1]])[c(7:30, 79:103)]
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
#sapply(seq_along(outDir), function(w) {
# system(paste0("[ -d ", outDir[w], " ] || mkdir ", outDir[w]))
#})
#sapply(seq_along(outDir), function(w) {
#  mclapply(seq_along(pop_name), function(x) {
#    system(paste0("[ -d ", outDir_list[[w]][x], " ] || mkdir ", outDir_list[[w]][x]))
#  }, mc.cores = length(pop_name), mc.preschedule = F)
#})
#sapply(seq_along(outDir), function(w) {
#  mclapply(seq_along(pop_name), function(x) {
#    system(paste0("[ -d ", plotDir_list[[w]][x], " ] || mkdir ", plotDir_list[[w]][x]))
#  }, mc.cores = length(pop_name), mc.preschedule = F)
#})

# For each population, divide features into quantiles based on decreasing orderingFactor
for(x in 1:length(featuresNLR_pop_list)) {
  print(pop_name[x])
  NLR_quants <- data.frame(featuresNLR_pop_list[[x]],
                           NLR_quantile = as.character(""),
                           stringsAsFactors = F)
  mclapply(seq_along(orderingFactor), function(w) {
    print(orderingFactor[w])
    # Assign 0s to NA values only for coverage data
    if(grepl("_in_", orderingFactor[w])) {
      NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])][which(is.na(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]))] <- 0
    }
    if(grepl("all", orderingFactor[w])) {
    # 2 quantiles for population genetics stats quantiles
      quantiles <- 2
    } else {
    # 4 quantiles for coverage and cM/Mb quantiles
     quantiles <- 4
    }
    quantilesStats <- data.frame()
    for(k in 1:quantiles) {
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
      if(k < quantiles) {
        NLR_quants[ !is.na(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) >  1-(k/quantiles), ]$NLR_quantile <- paste0("Quantile ", k)
      } else {
      # Final quantile should span 0 to, e.g., 0.25 proportions of features
        NLR_quants[ !is.na(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) <= 1-((k-1)/quantiles) &
                    percent_rank(NLR_quants[,which(colnames(NLR_quants) == orderingFactor[w])]) >= 1-(k/quantiles), ]$NLR_quantile <- paste0("Quantile ", k)
      }
      write.table(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),],
                  file = paste0(outDir_list[[w]][x],
                                "quantile", k, "_of_", quantiles,
                                "_by_", orderingFactor[w],
                                "_of_NLR_genes_in_",
                                paste0(substring(featureName, first = 10, last = 16),
                                       collapse = "_"), "_",
                                substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                  quote = FALSE, sep = "\t", row.names = FALSE)
      stats <- data.frame(quantile = as.integer(k),
                          n = as.integer(dim(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),])[1]),
                          mean_width = as.integer(round(mean(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),]$width, na.rm = T))),
                          total_width = as.integer(sum(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),]$width, na.rm = T)),
                          mean_orderingFactor = as.numeric(mean(NLR_quants[NLR_quants$NLR_quantile == paste0("Quantile ", k),][,which(colnames(NLR_quants) == orderingFactor[w])], na.rm = T)))
      quantilesStats <- rbind(quantilesStats, stats)
    }
    write.table(quantilesStats,
                file = paste0(outDir_list[[w]][x],
                              "summary_", quantiles, "quantiles",
                              "_by_", orderingFactor[w],
                              "_of_NLR_genes_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(NLR_quants,
                file = paste0(outDir_list[[w]][x],
                              "features_", quantiles, "quantiles",
                              "_by_", orderingFactor[w],
                              "_of_NLR_genes_in_",
                              paste0(substring(featureName, first = 10, last = 16),
                                     collapse = "_"), "_",
                              substring(featureName[1][1], first = 18), "_", pop_name[x], ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }, mc.cores = length(orderingFactor), mc.preschedule = F)
}
