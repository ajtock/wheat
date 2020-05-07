#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 100G -c 48 "/applications/R/R-3.3.2/bin/Rscript permTest_peaks_vs_others.R DMC1_Rep1_ChIP 10000 'euchromatin' 'A'"

args <- commandArgs(trailingOnly = TRUE)
libNameChIP <- args[1]
perms <- as.numeric(args[2])
region <- args[3]
genomeName <- args[4]

library(regioneR)

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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
genomeGR <- genomeGR[grep(genomeName,
                          seqnames(genomeGR))@values]

# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin or genomewide")
}

# Import peaks as GRanges object
load(paste0("../../../", libNameChIP,
            "_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
peaksGR <- peaksGR[grep(genomeName,
                        seqnames(peaksGR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = peaksGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
peaksGR <- peaksGR[-subjectHits(mask_peaks_overlap)]
strand(peaksGR) <- "*"
print("***********peaks***********")
print(peaksGR)

# ASY1
load("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/ASY1_CS_Rep1_ChIP_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
ASY1GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
ASY1GR <- ASY1GR[grep(genomeName,
                      seqnames(ASY1GR))@values]
# Subset to include only those not overlapping masked region (e.g., heterochromatin)
mask_peaks_overlap <- findOverlaps(query = maskGR,
                                   subject = ASY1GR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
ASY1GR <- ASY1GR[-subjectHits(mask_peaks_overlap)]
strand(ASY1GR) <- "*"
print("***********ASY1***********")
print(ASY1GR)

# H3K4me3
load("/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K4me3_Rep1_ChIP_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K4me3GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K4me3GR <- H3K4me3GR[grep(genomeName,
                            seqnames(H3K4me3GR))@values]
# Subset to include only those not overlapping masked region
mask_H3K4me3_overlap <- findOverlaps(query = maskGR,
                                     subject = H3K4me3GR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
H3K4me3GR <- H3K4me3GR[-subjectHits(mask_H3K4me3_overlap)]
strand(H3K4me3GR) <- "*"
print("***********H3K4me3***********")
print(H3K4me3GR)

# H3K9me2
load("/home/ajt200/analysis/wheat/H3K9me2/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K9me2_Rep1_ChIP_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K9me2GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K9me2GR <- H3K9me2GR[grep(genomeName,
                            seqnames(H3K9me2GR))@values]
# Subset to include only those not overlapping masked region
mask_H3K9me2_overlap <- findOverlaps(query = maskGR,
                                     subject = H3K9me2GR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
H3K9me2GR <- H3K9me2GR[-subjectHits(mask_H3K9me2_overlap)]
strand(H3K9me2GR) <- "*"
print("***********H3K9me2***********")
print(H3K9me2GR)

# H3K27me1
load("/home/ajt200/analysis/wheat/H3K27me1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K27me1_Rep1_ChIP_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K27me1GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K27me1GR <- H3K27me1GR[grep(genomeName,
                              seqnames(H3K27me1GR))@values]
# Subset to include only those not overlapping masked region
mask_H3K27me1_overlap <- findOverlaps(query = maskGR,
                                      subject = H3K27me1GR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
H3K27me1GR <- H3K27me1GR[-subjectHits(mask_H3K27me1_overlap)]
strand(H3K27me1GR) <- "*"
print("***********H3K27me1***********")
print(H3K27me1GR)

# H3K27me3
load("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K27me3/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K27me3_ChIP_SRR6350666_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K27me3GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K27me3GR <- H3K27me3GR[grep(genomeName,
                              seqnames(H3K27me3GR))@values]
# Subset to include only those not overlapping masked region
mask_H3K27me3_overlap <- findOverlaps(query = maskGR,
                                      subject = H3K27me3GR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
H3K27me3GR <- H3K27me3GR[-subjectHits(mask_H3K27me3_overlap)]
strand(H3K27me3GR) <- "*"
print("***********H3K27me3***********")
print(H3K27me3GR)

# H3K36me3
load("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K36me3/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K36me3_ChIP_SRR6350670_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K36me3GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K36me3GR <- H3K36me3GR[grep(genomeName,
                              seqnames(H3K36me3GR))@values]
# Subset to include only those not overlapping masked region
mask_H3K36me3_overlap <- findOverlaps(query = maskGR,
                                      subject = H3K36me3GR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
H3K36me3GR <- H3K36me3GR[-subjectHits(mask_H3K36me3_overlap)]
strand(H3K36me3GR) <- "*"
print("***********H3K36me3***********")
print(H3K36me3GR)

# H3K9ac
load("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/H3K9ac/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K9ac_ChIP_SRR6350667_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K9acGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K9acGR <- H3K9acGR[grep(genomeName,
                          seqnames(H3K9acGR))@values]
# Subset to include only those not overlapping masked region
mask_H3K9ac_overlap <- findOverlaps(query = maskGR,
                                    subject = H3K9acGR,
                                    type = "any",
                                    select = "all",
                                    ignore.strand = TRUE)
H3K9acGR <- H3K9acGR[-subjectHits(mask_H3K9ac_overlap)]
strand(H3K9acGR) <- "*"
print("***********H3K9ac***********")
print(H3K9acGR)

# H2AZ
load("/home/ajt200/analysis/wheat/H2AZ/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H2AZ_Rep1_ChIP_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H2AZGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H2AZGR <- H2AZGR[grep(genomeName,
                      seqnames(H2AZGR))@values]
# Subset to include only those not overlapping masked region
mask_H2AZ_overlap <- findOverlaps(query = maskGR,
                                  subject = H2AZGR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
H2AZGR <- H2AZGR[-subjectHits(mask_H2AZ_overlap)]
strand(H2AZGR) <- "*"
print("***********H2AZ***********")
print(H2AZGR)

# H3K4me1
load("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/H3K4me1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K4me1_Rep1_ChIP_SRR8126618_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K4me1GR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K4me1GR <- H3K4me1GR[grep(genomeName,
                            seqnames(H3K4me1GR))@values]
# Subset to include only those not overlapping masked region
mask_H3K4me1_overlap <- findOverlaps(query = maskGR,
                                     subject = H3K4me1GR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
H3K4me1GR <- H3K4me1GR[-subjectHits(mask_H3K4me1_overlap)]
strand(H3K4me1GR) <- "*"
print("***********H3K4me1***********")
print(H3K4me1GR)

# H3K27ac
load("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/H3K27ac/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/H3K27ac_Rep1_ChIP_SRR8126621_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K27acGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
H3K27acGR <- H3K27acGR[grep(genomeName,
                            seqnames(H3K27acGR))@values]
# Subset to include only those not overlapping masked region
mask_H3K27ac_overlap <- findOverlaps(query = maskGR,
                                     subject = H3K27acGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
H3K27acGR <- H3K27acGR[-subjectHits(mask_H3K27ac_overlap)]
strand(H3K27acGR) <- "*"
print("***********H3K27ac***********")
print(H3K27acGR)

# MNase
load("/home/ajt200/analysis/wheat/MNase/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/MNase_Rep1_rangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
MNaseGR <- rangerPeaksGRmergedOverlaps
rangerPeaksGRmergedOverlaps <- NULL
MNaseGR <- MNaseGR[grep(genomeName,
                        seqnames(MNaseGR))@values]
# Subset to include only those not overlapping masked region
mask_MNase_overlap <- findOverlaps(query = maskGR,
                                   subject = MNaseGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
MNaseGR <- MNaseGR[-subjectHits(mask_MNase_overlap)]
strand(MNaseGR) <- "*"
print("***********MNase***********")
print(MNaseGR)


# genes
genes <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA.gff3",
                    colClasses = c(NA,
                                   rep("NULL", 2),
                                   rep(NA, 2),
                                   "NULL", NA, "NULL", NA))
colnames(genes) <- c("chr", "start", "end", "strand", "geneID")
genes <- genes[genes$chr != "chrUn",]
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = genes$strand)
genesGR <- genesGR[grep(genomeName,
                        seqnames(genesGR))@values]
# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = maskGR,
                                   subject = genesGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
genesGR <- genesGR[-subjectHits(mask_genes_overlap)]
# Retain strand information until after obtaining promoters, etc.

# Obtain 1000-bp gene promoters
promotersGR <- promoters(genesGR, upstream = 1000, downstream = 0)
strand(promotersGR) <- "*"
print(promotersGR)

# Obtain regions immediately downstream of gene TSSs (TSS to TSS+499 bp)
TSSsGR <- promoters(genesGR, upstream = 0, downstream = 500)
strand(TSSsGR) <- "*"
print(TSSsGR)

# Obtain regions immediately upstream of gene TTSs (TTS to TTS-499 bp)
source("/projects/ajt200/Rfunctions/TTSplus.R")
TTSsGR <- TTSplus(genesGR, upstream = 499, downstream = 0)
strand(TTSsGR) <- "*"
print(TTSsGR)

# Obtain 1000-bp gene terminators
terminatorsGR <- TTSplus(genesGR, upstream = -1, downstream = 1000)
strand(terminatorsGR) <- "*"
print(terminatorsGR)

# Remove strand information from genesGR
strand(genesGR) <- "*"
print(genesGR)

# NLRs
NLRs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/",
                          "NLRs_Steuernagel_Wulff_2020_Plant_Physiol/NLR_genes_complete_representative_mRNA.gff3"),
                   colClasses = c(NA,
                                  rep("NULL", 2),
                                  rep(NA, 2),
                                  "NULL", NA, "NULL", NA))
NLRs$V9 <- sub(pattern = "\\.\\d+", replacement = "",
               x = NLRs$V9)
colnames(NLRs) <- c("chr", "start", "end", "strand", "geneID")
NLRs <- NLRs[NLRs$chr != "chrUn",]
geneIDs <- sub(pattern = "\\.\\d+", replacement = "",
               x = genes$geneID)
NLRs <- NLRs[NLRs$geneID %in% geneIDs,]
NLRsGR <- GRanges(seqnames = NLRs$chr,
                  ranges = IRanges(start = NLRs$start,
                                   end = NLRs$end),
                  strand = NLRs$strand)
NLRsGR <- NLRsGR[grep(genomeName,
                      seqnames(NLRsGR))@values]
# Subset to include only those not overlapping masked region
mask_NLRs_overlap <- findOverlaps(query = maskGR,
                                  subject = NLRsGR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
NLRsGR <- NLRsGR[-subjectHits(mask_NLRs_overlap)]
# Retain strand information until after obtaining promoters, etc.

# Obtain 1000-bp gene promoters
NLRpromotersGR <- promoters(NLRsGR, upstream = 1000, downstream = 0)
strand(NLRpromotersGR) <- "*"
print(NLRpromotersGR)

# Obtain regions immediately downstream of gene TSSs (TSS to TSS+499 bp)
NLRTSSsGR <- promoters(NLRsGR, upstream = 0, downstream = 500)
strand(NLRTSSsGR) <- "*"
print(NLRTSSsGR)

# Obtain regions immediately upstream of gene TTSs (TTS to TTS-499 bp)
source("/projects/ajt200/Rfunctions/TTSplus.R")
NLRTTSsGR <- TTSplus(NLRsGR, upstream = 499, downstream = 0)
strand(NLRTTSsGR) <- "*"
print(NLRTTSsGR)

# Obtain 1000-bp gene terminators
NLRterminatorsGR <- TTSplus(NLRsGR, upstream = -1, downstream = 1000)
strand(NLRterminatorsGR) <- "*"
print(NLRterminatorsGR)

# Remove strand information from NLRsGR
strand(NLRsGR) <- "*"
print(NLRsGR)

# meio
# Note: these two sets of meiotic genes share 271 common genes that are assigned to a chromosome
meio1 <- read.table("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/Table_S4_meiotic_GO_genes.tsv",
                    header = T, stringsAsFactors = F)$Gene.ID
meio2 <- read.table("/home/ajt200/analysis/wheat/RNAseq_meiocyte_Alabdullah_Moore_2019_FrontPlantSci/Table_S4_meiotic_gene_orthologs.tsv",
                    header = T, sep = "\t", stringsAsFactors = F)$Gene.ID
meio <- unique(c(meio1, meio2))
print(length(meio))
#[1] 1063
meio <- genes[sub(pattern = "\\.\\d+", replacement = "",
                  x = genes$geneID) %in% meio,]
print(dim(meio))
#[1] 1059    5
# Redundant but does no harm
meio <- meio[meio$chr != "chrUn",]
print(dim(meio))
#[1] 1059    5
# Convert into GRanges
meioGR <- GRanges(seqnames = meio$chr,
                  ranges = IRanges(start = meio$start,
                                   end = meio$end),
                  strand = meio$strand)
meioGR <- meioGR[grep(genomeName,
                      seqnames(meioGR))@values]
# Subset to include only those not overlapping masked region
mask_meio_overlap <- findOverlaps(query = maskGR,
                                  subject = meioGR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
meioGR <- meioGR[-subjectHits(mask_meio_overlap)]
# Retain strand information until after obtaining promoters, etc.

# Obtain 1000-bp gene promoters
meiopromotersGR <- promoters(meioGR, upstream = 1000, downstream = 0)
strand(meiopromotersGR) <- "*"
print(meiopromotersGR)

# Obtain regions immediately downstream of gene TSSs (TSS to TSS+499 bp)
meioTSSsGR <- promoters(meioGR, upstream = 0, downstream = 500)
strand(meioTSSsGR) <- "*"
print(meioTSSsGR)

# Obtain regions immediately upstream of gene TTSs (TTS to TTS-499 bp)
source("/projects/ajt200/Rfunctions/TTSplus.R")
meioTTSsGR <- TTSplus(meioGR, upstream = 499, downstream = 0)
strand(meioTTSsGR) <- "*"
print(meioTTSsGR)

# Obtain 1000-bp gene terminators
meioterminatorsGR <- TTSplus(meioGR, upstream = -1, downstream = 1000)
strand(meioterminatorsGR) <- "*"
print(meioterminatorsGR)

# Remove strand information from meioGR
strand(meioGR) <- "*"
print(meioGR)

## TEs
#TEdir <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/"
#TEs <- read.table(paste0(TEdir,
#                         "iwgsc_refseqv1.0_TransposableElements_2017Mar13.bed"),
#                         header = F)
#TEsGR <- GRanges(seqnames = TEs$V1,
#                 ranges = IRanges(start = TEs$V2+1,
#                                  end = TEs$V3),
#                 strand = "*")
#TEsGR <- TEsGR[grep(genomeName,
#                    seqnames(TEsGR))@values]
## Subset to include only those not overlapping masked region
#mask_TEs_overlap <- findOverlaps(query = maskGR,
#                                 subject = TEsGR,
#                                 type = "any",
#                                 select = "all",
#                                 ignore.strand = TRUE)
#TEsGR <- TEsGR[-subjectHits(mask_TEs_overlap)]
#print(TEsGR)
#
#otherNames2 <- c(
#                "CACTA_DTC",
#                "Harbinger_DTH",
#                "hAT_DTA",
#                "Helitrons_DHH",
#                "Mariner_DTT",
#                "Mutator_DTM",
#                "Unclassified_class_2_DXX",
#                "Unclassified_with_TIRs_DTX",
#                "MITE_DMI",
#                "Copia_LTR_RLC",
#                "Gypsy_LTR_RLG",
#                "LINE_RIX",
#                "SINE_SIX",
#                "Unclassified_LTR_RLX",
#                "Unclassified_repeats_XXX"
#               )
#TEfamNames <- c(
#                "CACTA_DTC",
#                "Harbinger_DTH",
#                "hAT_DTA",
#                "Helitrons_DHH",
#                "Mariner_DTT",
#                "Mutator_DTM",
#                "Unclassified_class_2_DXX",
#                "Unclassified_with_TIRs_DTX",
#                "MITE_DMI",
#                "Copia_LTR_RLC",
#                "Gypsy_LTR_RLG",
#                "LINE_RIX",
#                "SINE_SIX",
#                "Unclassified_LTR_RLX",
#                "Unclassified_repeats_XXX"
#               )
#
#TEdir2 <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/superfamilies/"
#
#### TEs
#othersGRL2 <- lapply(seq_along(TEfamNames), function(x) {
#  TEs <- read.table(paste0(TEdir,
#                           "iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
#                           TEfamNames[x], ".bed"), header = F)
#  TEsGRtmp <- GRanges(seqnames = TEs$V1,
#                      ranges = IRanges(start = TEs$V2+1,
#                                       end = TEs$V3),
#                      strand = "*")
#  TEsGRtmp <- TEsGRtmp[grep(genomeName,
#                            seqnames(TEsGRtmp))@values]
#  # Subset to include only those not overlapping masked region
#  mask_TEstmp_overlap <- findOverlaps(query = maskGR,
#                                      subject = TEsGRtmp,
#                                      type = "any",
#                                      select = "all",
#                                      ignore.strand = TRUE)
#  TEsGRtmp[-subjectHits(mask_TEstmp_overlap)]
#})

# Create vector of other-feature names
otherNames <- c(
                "ASY1",
                "H3K4me3",
                "H3K9me2",
                "H3K27me1",
                "H3K27me3",
                "H3K36me3",
                "H3K9ac",
                "H2AZ",
                "H3K4me1",
                "H3K27ac",
                "MNase",
                "genes",
                "promoters",
                "TSSsPlus500bp",
                "TTSsMinus500bp",
                "terminators",
                "NLRs",
                "NLRpromoters",
                "NLRTSSsPlus500bp",
                "NLRTTSsMinus500bp",
                "NLRterminators",
                "meio",
                "meiopromoters",
                "meioTSSsPlus500bp",
                "meioTTSsMinus500bp",
                "meioterminators"
               ) 
# Create GRangesList of other features
othersGRL <- c(
               "ASY1GR" = ASY1GR,
               "H3K4me3GR" = H3K4me3GR,
               "H3K9me2GR" = H3K9me2GR,
               "H3K27me1GR" = H3K27me1GR,
               "H3K27me3GR" = H3K27me3GR,
               "H3K36me3GR" = H3K36me3GR,
               "H3K9acGR" = H3K9acGR,
               "H2AZGR" = H2AZGR,
               "H3K4me1GR" = H3K4me1GR,
               "H3K27acGR" = H3K27acGR,
               "MNaseGR" = MNaseGR,
               "genesGR" = genesGR,
               "promotersGR" = promotersGR,
               "TSSsGR" = TSSsGR,
               "TTSsGR" = TTSsGR,
               "terminatorsGR" = terminatorsGR,
               "NLRsGR" = NLRsGR,
               "NLRpromotersGR" = NLRpromotersGR,
               "NLRTSSsGR" = NLRTSSsGR,
               "NLRTTSsGR" = NLRTTSsGR,
               "NLRterminatorsGR" = NLRterminatorsGR,
               "meioGR" = meioGR,
               "meiopromotersGR" = meiopromotersGR,
               "meioTSSsGR" = meioTSSsGR,
               "meioTTSsGR" = meioTTSsGR,
               "meioterminatorsGR" = meioterminatorsGR
              ) 

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(845934)
ptPeaksOtherPerChrom <- lapply(seq_along(othersGRL), function(x) {
  permTest(A = peaksGR,
           B = othersGRL[[x]],
           genome = genomeGR,
           mask = maskGR,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE,
           per.chromosome = TRUE,
           evaluate.function = numOverlaps,
           count.once = TRUE,
           ntimes = perms,
           mc.set.seed = FALSE,
           mc.cores = detectCores())
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0("permTest_", as.character(perms), "perms_",
                   libNameChIP, "_peaks_vs_others_in_",
                   genomeName, "genome_", region,
                   ".RData"))

# Summarise results in a table
featureName <- NULL
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksOtherPerChrom)) {
  featureNamei <- print(otherNames[i])
  featureName <- c(featureName, featureNamei)
  noOfFeaturesi <- print(length(othersGRL[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksOtherPerChromDataFrame <- cbind(featureName, noOfFeatures, expected, observed, pval, zscore)
colnames(ptPeaksOtherPerChromDataFrame) <- c("feature", "n", "expected", "observed", "pval", "zscore")
write.table(ptPeaksOtherPerChromDataFrame,
            file = paste0("permTest_", as.character(perms), "perms_",
                          libNameChIP, "_peaks_vs_others_in_",
                          genomeName, "genome_", region,
                          "_DataFrame.txt"),
            sep = "\t", quote = F, row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(paste0(plotDir, otherNames[i],
             "_permTest_", as.character(perms), "perms_",
             libNameChIP, "_peaks_in_",
             genomeName, "genome_", region,
             "_perChrom.pdf"),
      width = 10, height = 7)
  plot(ptPeaksOtherPerChrom[[i]], main = paste0(libNameChIP, " peaks vs ", otherNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = othersGRL[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = othersGRL[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = othersGRL[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(paste0(plotDir, otherNames[i],
             "_localZscore_permTest_", as.character(perms), "perms_",
             libNameChIP, "_peaks_in_",
             genomeName, "genome_", region,
             "_w1kb_s50bp_w10kb_s500bp_w",
             win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(libNameChIP, " peaks vs ", otherNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(libNameChIP, " peaks vs ", otherNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(libNameChIP, " peaks vs ", otherNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(libNameChIP, " peaks vs ", otherNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(libNameChIP, " peaks vs ", otherNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(libNameChIP, " peaks vs ", otherNames[i], " (~", win, "-kb shift)"))
  dev.off()
}
