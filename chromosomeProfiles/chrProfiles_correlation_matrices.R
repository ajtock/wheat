#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 28.04.2020

# Create and plot correlation matrices of chromosome profiles for
# log2(ChIP/input), MNase, DNase, input, RNA-seq, DNA methylation,
# and NLR-encoding, defense response and meiotic genes,
# genome-wide, subgenome-wide or separated into the following chromosome compartments a.k.a. partitions):
# 1. R1 and R3 (distal)
# 2. R2a and R2b (interstitial)
# 3. C (proximal)
# 4. heterochromatin (interstitial and proximal)
# 5. centromeres (defined by IWGSC (2018) Science 361 using CENH3 ChIP-seq data from Guo et al. (2016) PLOS Genet. 12)

# Usage:
# ./chrProfiles_correlation_matrices.R 1Mb 1000000 distal 'A,B,D' smoothed

#winName <- "1Mb"
#winSize <- 1000000
#region <- "distal"
#genomeName <- unlist(strsplit("A,B,D",
#                              split = ","))
#smoothing <- "smoothed"

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])
region <- args[3]
genomeName <- unlist(strsplit(args[4],
                              split = ","))
smoothing <- args[5]

plotDir <- paste0("plots/Figure1/correlation_matrices/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

library(GenomicRanges)
library(Hmisc) # includes rcorr() function which computes significance levels for Pearson and Spearman correlations
library(reshape)
library(ggplot2)
library(ggcorrplot)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                               "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,2])
centromereEnd <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                             "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = chrStart,
                                     end = chrLens),
                    strand = "*")
if(length(genomeName) == 1) {
  genomeGR <- genomeGR[grep(genomeName,
                            seqnames(genomeGR))@values]
}

# Define region to be analysed
if(region == "distal") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "interstitial") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(chrPartitions$R1_R2a+1,
                                                 chrPartitions$C_R2b),
                                       end = c(chrPartitions$R2a_C,
                                               chrPartitions$R2b_R3-1)),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "proximal") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  if(length(genomeName) == 1) {
    regionGR <- regionGR[grep(genomeName,
                              seqnames(regionGR))@values]
  }
} else {
  stop("region is not distal, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "distal") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "interstitial") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 3),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2a_C+1,
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrPartitions$C_R2b-1,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "proximal") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else if(region == "genomewide") {
  maskGR <- GRanges()
  if(length(genomeName) == 1) {
    maskGR <- maskGR[grep(genomeName,
                          seqnames(maskGR))@values]
  }
} else {
  stop("region is not distal, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

ChIPDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/log2ChIPcontrol/"
otherDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/MNase_DNase_input_RNAseq/"
cMMbDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/cMMb/"
geneDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/genes/"
DNAmethDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/DNAmeth/"

paths <- c(
           paste0(ChIPDir, "log2_DMC1_Rep1_ChIP_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_ASY1_CS_Rep1_ChIP_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(cMMbDir, "cMMb_iwgsc_refseqv1.0_mapping_data_minInterMarkerDist200bp_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K4me1_Rep1_ChIP_SRR8126618_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K4me3_Rep1_ChIP_MNase_Rep1_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K9ac_ChIP_SRR6350667_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K27ac_Rep1_ChIP_SRR8126621_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K27me3_ChIP_SRR6350666_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K36me3_ChIP_SRR6350670_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H2AZ_Rep1_ChIP_MNase_Rep1_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K9me2_Rep1_ChIP_MNase_Rep1_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_H3K27me1_Rep1_ChIP_MNase_Rep1_per_", winName, "_", smoothing, ".txt"),
           paste0(DNAmethDir, "BSseq_Rep8a_SRR6792678_per_", winName, "_", smoothing, ".txt"),
           paste0(ChIPDir, "log2_CENH3_ChIP_SRR1686799_H3_input_SRR6350669_per_", winName, "_", smoothing, ".txt"),
           paste0(otherDir, "MNase_Rep1_per_", winName, "_", smoothing, ".txt"),
           paste0(otherDir, "WT_RNAseq_Rep1_ERR2402974_per_", winName, "_", smoothing, ".txt"),
           paste0(geneDir, "gene_frequency_per_", winName, "_", smoothing, ".txt"),
           paste0(geneDir, "defense_response_gene_frequency_per_", winName, "_", smoothing, ".txt"),
           paste0(geneDir, "NLR_gene_frequency_per_", winName, "_", smoothing, ".txt"),
           paste0(geneDir, "meiotic_gene_frequency_per_", winName, "_", smoothing, ".txt")
)

profileNames <- c(
                  "DMC1",
                  "ASY1",
                  "cM/Mb",
                  "H3K4me1",
                  "H3K4me3",
                  "H3K9ac",
                  "H3K27ac",
                  "H3K27me3",
                  "H3K36me3",
                  "H2A.Z",
                  "H3K9me2",
                  "H3K27me1",
                  "mCG",
                  "mCHG",
                  "mCHH",
                  "CENH3",
                  "MNase",
                  "RNA-seq (meiocytes)",
                  "Genes",
                  "Defense response genes",
                  "NLR-encoding genes",
                  "Meiotic genes"
)

profiles <- lapply(seq_along(paths), function(x) {
  read.table(paths[x], header = T)
})
profilesNew <- lapply(seq_along(profiles), function(x) {
  # Make columns consistent across log2ChIPinput,
  # cMMb, MNase, DNase, RNA-seq and gene data sets
  if(dim(profiles[[x]])[2] >= 3 &
     dim(profiles[[x]])[2] < 6) {
    data.frame(chr = profiles[[x]][,1],
               window = profiles[[x]][,2],
               value = profiles[[x]][,dim(profiles[[x]])[2]],
               stringsAsFactors = F)
  # Separate DNA methylation contexts and make columns consistent as above
  } else if(dim(profiles[[x]])[2] == 6) {
    DNAmethList <- list()
#    for(y in 3:dim(profiles[[x]])[2]) {
    # Exlcude average over all 3 contexts
    for(y in 3:5) {
      print(y)
      DNAmethList[[y]] <- data.frame(chr = profiles[[x]][,1],
                                     window = profiles[[x]][,2],
                                     value = profiles[[x]][,y],
                                     stringsAsFactors = F)
    }
    # Remove empty ("NULL") list elements
    DNAmethList[-1:-2]
  }
})

# DNAmethList is one list of 3 elements within the 20-element list profilesNew
# To make each of these 3 elements its own element within a 22-element list profilesNew2:
profilesNew2 <- profilesNew
for(x in seq_along(profilesNew)) {
  if(class(profilesNew[[x]]) == "list") {
    listLength <- length(profilesNew[[x]])
    for(z in (length(profilesNew) + listLength - 1) : (x+listLength) ) {
      profilesNew2[[z]] <- profilesNew[[z - listLength + 1]]
    }
    for(y in 1:3) {
      profilesNew2[[x + y - 1]] <- profilesNew[[x]][[y]]
    }
  }
}

profiles <- profilesNew2 

profilesGR <- lapply(seq_along(profiles), function(x) {
  GRanges(seqnames = profiles[[x]]$chr,
          ranges = IRanges(start = profiles[[x]]$window,
                           end = profiles[[x]]$window+winSize-1),
          strand = "*",
          value = profiles[[x]]$value)
}) 
# Redefine end coordinate of last window in each chromosome
# to be equal to the length of that chromosome
profilesGR_chrs <- lapply(seq_along(profilesGR), function(x) {
  profilesGRx <- lapply(seq_along(chrs), function(i) {
    profilesGRx_chr <- profilesGR[[x]][seqnames(profilesGR[[x]]) == chrs[i]]
    end(profilesGRx_chr)[length(profilesGRx_chr)] <- chrLens[i]
    profilesGRx_chr
  })
  profilesGRx
})

# Combine GRanges list elements (1 element for each chromosome)
# in one GRanges object for each profile
profilesGR <- lapply(seq_along(profilesGR_chrs), function(x) {
  do.call(c, profilesGR_chrs[[x]])
})

# Subset to include only those windows within the genomeName-subgenome
if(length(genomeName) == 1) {
  profilesGR <- lapply(seq_along(profilesGR), function(x) {
    profilesGR[[x]][grep(genomeName,
                         seqnames(profilesGR[[x]]))@values]
  })
}
# Subset to include only those windows not overlapping masked region (e.g., heterochromatin)
profilesGR <- lapply(seq_along(profilesGR), function(x) {
  mask_profilesGRx_overlap <- findOverlaps(query = maskGR,
                                           subject = profilesGR[[x]],
                                           type = "any",
                                           select = "all",
                                           ignore.strand = TRUE)
  profilesGR[[x]][-subjectHits(mask_profilesGRx_overlap)]
})

# Combine profiles into one data.frame in which each profile is a column
profilesVal <- lapply(seq_along(profilesGR), function(x) {
  profilesGR[[x]]$value
})

profilesDF <- as.data.frame(do.call(cbind, profilesVal),
                            stringsAsFactors = F)
colnames(profilesDF) <- profileNames

# Create correlation matrix
corMat <- round(cor(profilesDF,
                    method = "spearman",
                    use = "pairwise.complete.obs"),
                digits = 2)
# Create matrix of P-values for correlation matrix
corMatSig <- rcorr(as.matrix(profilesDF), type = "spearman")$P
corMatSig <- corMatSig[,-1]

# Set duplicates to NA
for(x in 1:dim(corMat)[1]) {
  corMat[x, x] <- NA
  if(x > 1) {
    corMat[x, 1:x-1] <- NA
  }
}
corMat <- corMat[,-1]

# Convert into reshape::melt formatted data.frame
# and remove duplicate pairs
corDat <- melt(corMat)
corDat <- corDat[-which(is.na(corDat[,3])),]

# Order the data.frame for plotting
profileNamesList <- as.list(profileNames)
names(profileNamesList) <- profileNames
levels(corDat$X1) <- rev(profileNamesList) 
levels(corDat$X2) <- profileNamesList[-1]

# Plot
ggObj <- ggplot(data = corDat,
                mapping = aes(X2, X1, fill = value)) +
  geom_tile() +
  geom_text(mapping = aes(X2, X1, label = value)) +
  scale_fill_gradient2(name = bquote("Spearman's" ~ italic(r[s])),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, breaks = seq(-1, 1, by = 0.4), limits = c(-1, 1)) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(barwidth = 14, barheight = 2,
                                title.position = "top", title.hjust = 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0, size = 14, colour = "black"),
        axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.justification = c(1, 0),
        legend.position = c(0.45, 0.1),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(5.5, 70.5, 5.5, 5.5), "pt"),
        plot.title = element_text(hjust = 0.5, size = 20, colour = "black")) +
  ggtitle(bquote(.(winSize/1e6) * "-Mb Spearman's" ~ italic(r[s]) ~ "for" ~
          .(paste0(genomeName, collapse = "-, ")) * "-genome" ~
          .(region) ~ "regions (" * .(smoothing) * ")"))
ggsave(paste0(plotDir,
              "Spearman_correlation_matrix_", winName,
              "_log2ChIPcontrol_cMMb_MNase_DNase_input_RNAseq_DNAmeth_genes_in_",
              paste0(genomeName, collapse = "_"), "_genome_", region, "_", smoothing, ".pdf"),
       plot = ggObj, height = 12, width = 12)
