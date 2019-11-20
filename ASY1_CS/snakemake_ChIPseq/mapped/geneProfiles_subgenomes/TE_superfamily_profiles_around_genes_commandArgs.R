#!/applications/R/R-3.4.0/bin/Rscript

# Profile TE density around genes and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./TE_superfamily_profiles_around_genes_commandArgs.R genes_in_Agenome genomewide 3500 2000 2kb 20

#featureName <- "genes_in_Agenome"
#region <- "genomewide"
#bodyLength <- 3500
#flankSize <- 2000
#flankName <- "2kb"
#winSize <- 20

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
region <- args[2]
bodyLength <- as.numeric(args[3])
flankSize <- as.numeric(args[4])
flankName <- as.character(args[5])
winSize <- as.numeric(args[6])

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)
library(parallel)

matDir <- paste0("matrices/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

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

# Load genes in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
genes <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                           substring(featureName, first = 10, last = 16), "_", region, ".bed"),
                    header = F)
colnames(genes) <- c("chr", "start", "end", "name", "score", "strand")
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start+1,
                                    end = genes$end),
                   strand = genes$strand,
                   number = genes$name)
genesGR <- genesGR[seqnames(genesGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = maskGR,
                                   subject = genesGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
if(length(mask_genes_overlap) > 0) {
  genesGR <- genesGR[-subjectHits(mask_genes_overlap)]
}
# Load ranLoc in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLoc <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                            substring(featureName, first = 10, last = 16), "_", region, "_randomLoci.bed"),
                     header = F)
colnames(ranLoc) <- c("chr", "start", "end", "name", "score", "strand")
ranLocGR <- GRanges(seqnames = ranLoc$chr,
                    ranges = IRanges(start = ranLoc$start+1,
                                     end = ranLoc$end),
                    strand = ranLoc$strand,
                    number = ranLoc$name)
ranLocGR <- ranLocGR[seqnames(ranLocGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_ranLoc_overlap <- findOverlaps(query = maskGR,
                                    subject = ranLocGR,
                                    type = "any",
                                    select = "all",
                                    ignore.strand = TRUE)
if(length(mask_ranLoc_overlap) > 0) {
  ranLocGR <- ranLocGR[-subjectHits(mask_ranLoc_overlap)]
}

# Load TE superfamily BED files
inDirSuperfams <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/superfamilies/"
superfamCode <- c("RLG",
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
superfamName <- c("Gypsy_LTR",
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

superfamListGR <- mclapply(seq_along(superfamName), function(x) {
  superfam <- read.table(paste0(inDirSuperfams,
                                "iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
                                superfamName[x], "_", superfamCode[x], ".bed"),
                         header = F)
  colnames(superfam) <- c("chr", "start", "end", "name", "score", "strand")
  superfamGR <- GRanges(seqnames = superfam$chr,
                        ranges = IRanges(start = superfam$start+1,
                                         end = superfam$end),
                        strand = "*",
                        number = superfam$name,
                        coverage = rep(1, dim(superfam)[1]))
  superfamGR <- superfamGR[seqnames(superfamGR) != "chrUn"]
  # Subset to include only those not overlapping masked region
  mask_superfam_overlap <- findOverlaps(query = maskGR,
                                        subject = superfamGR,
                                        type = "any",
                                        select = "all",
                                        ignore.strand = TRUE)
  if(length(mask_superfam_overlap) > 0) {
    superfamGR <- superfamGR[-subjectHits(mask_superfam_overlap)]
  }
  superfamGR
}, mc.cores = length(superfamName))

# Define matrix and column mean outfiles
outDF <- lapply(seq_along(superfamName), function(x) {
  list(paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_", featureName, "_", region,
              "_matrix_bin", winSize, "bp_flank", flankName, ".tab"),
       paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_", featureName, "_", region,
              "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, ".tab"))
})
outDFcolMeans <- lapply(seq_along(superfamName), function(x) {
  list(paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_", featureName, "_", region,
              "_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"),
       paste0(matDir, superfamName[x], "_", superfamCode[x],
              "_around_", featureName, "_", region,
              "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"))
})

# Run covMatrix() function on each feature GRanges object to obtain matrices
# containing normalised feature density values around target and random loci
mclapply(seq_along(superfamName), function(x) {
  covMatrix(signal = superfamListGR[[x]],
            feature = genesGR,
            ranLoc = ranLocGR,
            featureSize = bodyLength,
            flankSize = flankSize,
            winSize = winSize,
            outDF = outDF[[x]],
            outDFcolMeans = outDFcolMeans[[x]])
  print(paste0(superfamName[x], "_", superfamCode[x],
               "_around_", featureName, "_", region,
               " profile calculation complete"))
}, mc.cores = length(superfamName))
