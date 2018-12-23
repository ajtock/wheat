#!/applications/R/R-3.5.0/bin/Rscript

library(rtracklayer)
library(parallel)
library(regioneR)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrStart <- rep(1, times = length(chrs))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/"
outDir <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/"
outDirSuperfams <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/superfamiliesTest/"
outDirSubfams <- "/home/ajt200/analysis/wheat/featureProfiles/TEs/subfamilies/"

TEsGFF <- readGFF(paste0(inDir, "iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3"))

superfamCode <- c("RLG",
                  "RLC",
                  "RLX",
                  "RIX",
                  "SIX",
                  "DTC",
                  "DTM",
                  "DTX",
                  "DTH",
                  "DTT",
                  "DXX",
                  "DTA",
                  "DHH",
                  "XXX")
superfamName <- c("Gyspy_LTR",
                  "Copia_LTR",
                  "Unclassified_LTR",
                  "LINE",
                  "SINE",
                  "CACTA",
                  "Mutator",
                  "Unclassified_with_TIRs",
                  "Harbinger",
                  "Mariner",
                  "Unclassified_class_2",
                  "hAT",
                  "Helitrons",
                  "Unclassified_repeats")

# Create a list of TE superfamily DataFrames
superfamDFList <- mclapply(seq_along(superfamCode), function(x)
  {
    TEsGFF[grep(superfamCode[x], TEsGFF$compo),]
  }, mc.cores = length(superfamCode))

# Create a list of TE superfamilies, wherein each list element
# is a vector of TE subfamilies
subfamList <- mclapply(seq_along(superfamDFList), function(x)
  {
    # Extract the first subfamily specified in the
    # composition ("compo") column
    unique(gsub(" .*$", "", superfamDFList[[x]]$compo))
  }, mc.cores = length(superfamDFList))

# Convert superfamily DFs to BED and save
superfamBEDList <- mclapply(seq_along(superfamDFList), function(x)
  {
    data.frame(chr = as.character(superfamDFList[[x]]$seqid),
               start = as.integer(superfamDFList[[x]]$start-1),
               end = as.integer(superfamDFList[[x]]$end),
               name = as.character(seq(from = 1,
                                       to = dim(superfamDFList[[x]])[1],
                                       by = 1)),
               score = as.numeric(superfamDFList[[x]]$score),
               strand = as.character(superfamDFList[[x]]$strand))
  }, mc.cores = length(superfamDFList))

mclapply(seq_along(superfamBEDList), function(x)
  {
    write.table(superfamBEDList[[x]],
                file = paste0(outDirSuperfams,
                              "iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
                              superfamName[x], "_", superfamCode[x], ".bed"),
                quote = F, sep = "\t", row.names = F, col.names = F)
  }, mc.cores = length(superfamBEDList))

allBED <- data.frame(chr = as.character(TEsGFF$seqid),
                     start = as.integer(TEsGFF$start-1),
                     end = as.integer(TEsGFF$end),
                     name = as.character(seq(from = 1,
                                             to = dim(TEsGFF)[1],
                                             by = 1)),
                     score = as.numeric(TEsGFF$score),
                     strand = as.character(TEsGFF$strand))
write.table(allBED,
            file = paste0(outDir,
                          "iwgsc_refseqv1.0_TransposableElements_2017Mar13.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Convert superfamily DFs to GRanges
superfamGRList <- mclapply(seq_along(superfamDFList), function(x)
  {
    GRanges(seqnames = superfamDFList[[x]]$seqid,
            ranges = IRanges(start = superfamDFList[[x]]$start,
                             end = superfamDFList[[x]]$end),
            strand = superfamDFList[[x]]$strand)
  }, mc.cores = length(superfamDFList))

allGR <- GRanges(seqnames = TEsGFF$seqid,
                 ranges = IRanges(start = TEsGFF$start,
                                  end = TEsGFF$end),
                 strand = TEsGFF$strand)

# Generate coordinates for the same number of random loci
superfamRanLocGRList <- mclapply(seq_along(superfamGRList), function(x)
  {
    randomizeRegions(superfamGRList[[x]],
                     genome = genome,
                     per.chromosome = TRUE,
                     allow.overlaps = TRUE)
  }, mc.cores = length(superfamGRList))

allRanLocGR <- randomizeRegions(allGR,
                                genome = genome,
                                per.chromosome = TRUE,
                                allow.overlaps = TRUE)
                     
# Convert random loci GRanges to BED and save
superfamRanLocBEDList <- mclapply(seq_along(superfamRanLocGRList), function(x)
  {
    data.frame(chr = as.character(seqnames(superfamRanLocGRList[[x]])),
               start = as.integer(start(superfamRanLocGRList[[x]])-1),
               end = as.integer(end(superfamRanLocGRList[[x]])))
  }, mc.cores = length(superfamRanLocGRList))
mclapply(seq_along(superfamRanLocBEDList), function(x)
  {
    write.table(superfamRanLocBEDList[[x]],
                file = paste0(outDirSuperfams,
                              "iwgsc_refseqv1.0_TransposableElements_2017Mar13_superfamily_",
                              superfamName[x], "_", superfamCode[x], "_randomLoci.bed"),
                quote = F, sep = "\t", row.names = F, col.names = F)
  }, mc.cores = length(superfamRanLocBEDList))

allRanLocBED <- data.frame(chr = as.character(seqnames(allRanLocGR)),
                           start = as.integer(start(allRanLocGR)-1),
                           end = as.integer(end(allRanLocGR)))
write.table(allRanLocBED,
            file = paste0(outDir,
                          "iwgsc_refseqv1.0_TransposableElements_2017Mar13_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
