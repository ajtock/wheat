#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./stress_response_gene_frequency_chrProfiles.R 1Mb 1000000

#winName <- "1Mb"
#winSize <- 1000000

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])

library(rtracklayer)
library(parallel)
library(GenomicRanges)

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
outDir <- "/home/ajt200/analysis/wheat/chromosomeProfiles/genes/"
genesGFF <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"))
genesGFF$group <- sub(pattern = "\\.\\d+", replacement = "",
                      x = genesGFF$group)

# Load functional annotation in order to extract defensins
anno <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned.TAB",
                   header = T, quote = "\"", sep = "\t", check.names = F)
# Retain columns containing gene IDs and human-readable descriptions
anno <- data.frame(geneID = anno$`Gene-ID`,
                   description = anno$`Human-Readable-Description`,
                   GO = anno$`GO-IDs-(Description)-via-Interpro`,
                   stringsAsFactors = F)
# Replace gene model ID decimal suffix (e.g., ".1")
anno$geneID <- sub(pattern = "\\.\\d+", replacement = "",
                   x = anno$geneID)
# Replace "1G" with "2G" in gene IDs for consistency with v1.1
anno$geneID <- sub(pattern = "1G", replacement = "2G",
                   x = anno$geneID)

response_to_stress_indices <- which(grepl(pattern = "response to stress", x = anno$GO,
                              ignore.case = T))
defen_indices <- which(grepl(pattern = "defen", x = anno$GO,
                             ignore.case = T))
wounding_indices <- which(grepl(pattern = "wounding", x = anno$GO,
                                ignore.case = T))
response_to_heat_indices <- which(grepl(pattern = "response to heat", x = anno$GO,
                                        ignore.case = T))
response_to_water_indices <- which(grepl(pattern = "response to water", x = anno$GO,
                                         ignore.case = T))
desiccation_indices <- which(grepl(pattern = "desiccation", x = anno$GO,
                                   ignore.case = T))
indices <- sort(unique(c(response_to_stress_indices, defen_indices, wounding_indices,
                         response_to_heat_indices, response_to_water_indices, desiccation_indices))) 

# Retain only stress response genes
anno <- anno[indices,]

# Retain genes whose IDs are present in anno
genesGFF <- genesGFF[genesGFF$group %in% anno$geneID,]

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
#chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
#chrLens <- chrLens[-length(chrLens)]

# Define windows as GRanges object
windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  seqWindows <- seq(1, chrLens[i], by = winSize)
  windowsIR <- IRanges(start = seqWindows,
                       width = winSize)
  windowsIR <- windowsIR[-length(windowsIR)]
  windowsIR <- append(windowsIR,
                      IRanges(start = seqWindows[length(seqWindows)],
                              end = chrLens[i]))
  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = windowsIR,
                          strand = "*")
  print(chrWindowsGR)
  windowsGR <- append(windowsGR, chrWindowsGR)
}

featureProfile <- NULL
for(i in 1:length(chrs)) {
  # Count genes within windows
  chrgenes <- genesGFF[genesGFF$seqid == chrs[i],]
  chrgenesGR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chrgenes$start,
                                         end = chrgenes$end),
                        strand = chrgenes$strand)
  chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
  wingenes <- countOverlaps(chrWindowsGR,
                            chrgenesGR,
                            ignore.strand = T)
  chrProfile <- data.frame(chr = as.character(chrs[i]),
                           window = as.integer(start(chrWindowsGR)),
                           wingenes = as.integer(wingenes))
  featureProfile <- rbind(featureProfile, chrProfile)
}
write.table(featureProfile,
            file = paste0(outDir,
                          "stress_response_gene_frequency_per_", winName,
                          ".txt"))
