#!/applications/R/R-3.5.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 08.06.2020

# Profile SNP frequency around genes and random loci

# Wheat subgenome compartments (a.k.a. partitions):
# 1. R1 and R3 (euchromatin)
# 2. R2a and R2b (interstitial)
# 3. C (proximal)
# 4. heterochromatin (interstitial and proximal)
# 5. centromeres (defined by IWGSC (2018) Science 361 using CENH3 ChIP-seq data from Guo et al. (2016) PLOS Genet. 12)

# Usage via Condor submission system on node7:
# csmit -m 100G -c 13 "/applications/R/R-3.5.0/bin/Rscript ./15cvs_v_CS_WGS_SNP_profiles_around_genes_commandArgs.R genes_in_Agenome genomewide 3500 2000 2kb 20"

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

matDir <- paste0("matrices/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

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
genomeGR <- genomeGR[grep(substr(featureName, 10, 10),
                          seqnames(genomeGR))@values]

# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "interstitial") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(chrPartitions$R1_R2a+1,
                                                 chrPartitions$C_R2b),
                                       end = c(chrPartitions$R2a_C,
                                               chrPartitions$R2b_R3-1)),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "proximal") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else {
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "interstitial") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 3),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2a_C+1,
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrPartitions$C_R2b-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "proximal") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
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

# Load SNPs
SNPsDir <- "/home/ajt200/analysis/wheat/Brinton_Uauy_15cvs_v_CS_WGS_SNPs/all_SNPs_CS_position/snp_positions/"
SNPsList <- mclapply(seq_along(chrs), function(x) {
  read.table(paste0(SNPsDir,
                    "snp_positions_CS_ref_", chrs[x], ".tsv"),
             header = T)
}, mc.cores = length(chrs), mc.preschedule = F)
SNPs <- do.call(rbind, SNPsList)
colnames(SNPs) <- c("chr", "pos")
#SNPs$chr <- sub(pattern = "^", replacement = "chr", SNPs$chr)
## all SNPs
SNPsGR <- GRanges(seqnames = paste0("chr", SNPs$chr),
                  ranges = IRanges(start = SNPs$pos,
                                   end = SNPs$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs)[1]))

# Define matrix and column mean frequency outfile (mean profiles)
outDF <- list(paste0(matDir,
                     "15cvs_v_CS_WGS_all_SNPs_around_", featureName, "_", region,
                     "_matrix_bin", winSize, "bp_flank", flankName, ".tab"),
              paste0(matDir,
                     "15cvs_v_CS_WGS_all_SNPs_around_", featureName, "_", region,
                     "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, ".tab"))
outDFcolMeans <- list(paste0(matDir,
                             "15cvs_v_CS_WGS_all_SNPs_around_", featureName, "_", region,
                             "_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"),
                      paste0(matDir,
                             "15cvs_v_CS_WGS_all_SNPs_around_", featureName, "_", region,
                             "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing frequency values around target and random loci
covMatrix(signal = SNPsGR,
          feature = genesGR,
          ranLoc = ranLocGR,
          featureSize = bodyLength,
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0("SNP frequency around ", featureName, " ", region,
             " profile calculation complete"))

# Load SNPs and indels
SNPs_indelsDir <- "/home/ajt200/analysis/wheat/Brinton_Uauy_15cvs_v_CS_WGS_SNPs/all_SNPs_CS_position/snp_positions_indels/"
SNPs_indelsList <- mclapply(seq_along(chrs), function(x) {
  read.table(paste0(SNPs_indelsDir,
                    "snp_and_indel_positions_CS_ref_", chrs[x], ".tsv"),
             header = T)
}, mc.cores = length(chrs), mc.preschedule = F)
SNPs_indels <- do.call(rbind, SNPs_indelsList)
colnames(SNPs_indels) <- c("chr", "pos")
#SNPs_indels$chr <- sub(pattern = "^", replacement = "chr", SNPs_indels$chr)
## all SNPs_indels
SNPs_indelsGR <- GRanges(seqnames = paste0("chr", SNPs_indels$chr),
                         ranges = IRanges(start = SNPs_indels$pos,
                                          end = SNPs_indels$pos),
                         strand = "*",
                         coverage = rep(1, dim(SNPs_indels)[1]))

# Define matrix and column mean frequency outfile (mean profiles)
outDF <- list(paste0(matDir,
                     "15cvs_v_CS_WGS_all_SNPs_indels_around_", featureName, "_", region,
                     "_matrix_bin", winSize, "bp_flank", flankName, ".tab"),
              paste0(matDir,
                     "15cvs_v_CS_WGS_all_SNPs_indels_around_", featureName, "_", region,
                     "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, ".tab"))
outDFcolMeans <- list(paste0(matDir,
                             "15cvs_v_CS_WGS_all_SNPs_indels_around_", featureName, "_", region,
                             "_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"),
                      paste0(matDir,
                             "15cvs_v_CS_WGS_all_SNPs_indels_around_", featureName, "_", region,
                             "_ranLoc_matrix_bin", winSize, "bp_flank", flankName, "_colMeans.tab"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing frequency values around target and random loci
covMatrix(signal = SNPs_indelsGR,
          feature = genesGR,
          ranLoc = ranLocGR,
          featureSize = bodyLength,
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0("SNP and indel frequency around ", featureName, " ", region,
             " profile calculation complete"))
