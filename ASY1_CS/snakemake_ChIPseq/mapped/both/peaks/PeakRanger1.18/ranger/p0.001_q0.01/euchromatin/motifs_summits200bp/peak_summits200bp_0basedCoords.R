# Define peak coordinates as peak summit-100 bp to peak summit+100 bp
# so that all peaks have a common width of 201 bp for motif analysis
# Generate peak bed files (0-based start coordinates)

source("/projects/ajt200/Rfunctions/locus_midpoint.R")
library(GenomicRanges)
library(regioneR)

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
seqlevels(genome) <- sub("Chr", "", seqlevels(genome))

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/motifs_summits200bp/"

load(paste0(inDir,
            "REC8_HA_Rep1_armrangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_noMinWidth.RData"))
seqlevels(armrangerPeaksGR) <- sub("Chr", "", seqlevels(armrangerPeaksGR))
armrangerPeaksGR_summits <- GRanges(seqnames = seqnames(armrangerPeaksGR),
                                    ranges = IRanges(start = start(armrangerPeaksGR)+armrangerPeaksGR$summit0based,
                                                     end = start(armrangerPeaksGR)+armrangerPeaksGR$summit0based),
                                    strand = "*")
armrangerPeaksGR_summits200bp <- locMidpointFlank(armrangerPeaksGR_summits,
                                                  leftFlank = 100,
                                                  rightFlank = 100)
armrangerPeaks_summits200bp_0based <- cbind(as.numeric(seqnames(armrangerPeaksGR_summits200bp)),
                                            start(armrangerPeaksGR_summits200bp)-1,
                                            end(armrangerPeaksGR_summits200bp))
colnames(armrangerPeaks_summits200bp_0based) <- c("chr", "start", "end")
write.table(armrangerPeaks_summits200bp_0based,
            file = paste0(outDir, "armrangerPeaks_summits200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)

load(paste0(inDir,
            "REC8_HA_Rep1_perirangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_noMinWidth.RData"))
seqlevels(perirangerPeaksGR) <- sub("Chr", "", seqlevels(perirangerPeaksGR))
perirangerPeaksGR_summits <- GRanges(seqnames = seqnames(perirangerPeaksGR),
                                    ranges = IRanges(start = start(perirangerPeaksGR)+perirangerPeaksGR$summit0based,
                                                     end = start(perirangerPeaksGR)+perirangerPeaksGR$summit0based),
                                    strand = "*")
perirangerPeaksGR_summits200bp <- locMidpointFlank(perirangerPeaksGR_summits,
                                                   leftFlank = 100,
                                                   rightFlank = 100)
perirangerPeaks_summits200bp_0based <- cbind(as.numeric(seqnames(perirangerPeaksGR_summits200bp)),
                                             start(perirangerPeaksGR_summits200bp)-1,
                                             end(perirangerPeaksGR_summits200bp))
colnames(perirangerPeaks_summits200bp_0based) <- c("chr", "start", "end")
write.table(perirangerPeaks_summits200bp_0based,
            file = paste0(outDir, "perirangerPeaks_summits200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)


# Generate matched random loci

## arm
# mask pericentromeric regions
maskPeri <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))
seqlevels(maskPeri) <- sub("Chr", "", seqlevels(maskPeri))

armranLocGR_200bp <- randomizeRegions(armrangerPeaksGR_summits200bp,
                                      genome = genome,
                                      mask = c(maskPeri, armrangerPeaksGR_summits200bp),
                                      per.chromosome = T, allow.overlaps = T)
armranLoc_200bp_0based <- cbind(as.numeric(seqnames(armranLocGR_200bp)),
                                start(armranLocGR_200bp)-1, end(armranLocGR_200bp))
colnames(armranLoc_200bp_0based) <- c("chr", "start", "end")
write.table(armranLoc_200bp_0based, file = paste0(outDir, "armranLoc_200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)

## peri
# mask arm regions
maskArms <- toGRanges(data.frame(rep(chrs, 2), c(chrStart, pericenEnd), c(pericenStart, chrLens)))
seqlevels(maskArms) <- sub("Chr", "", seqlevels(maskArms))

periranLocGR_200bp <- randomizeRegions(perirangerPeaksGR_summits200bp,
                                       genome = genome,
                                       mask = c(maskArms, perirangerPeaksGR_summits200bp),
                                       per.chromosome = T, allow.overlaps = T)
periranLoc_200bp_0based <- cbind(as.numeric(seqnames(periranLocGR_200bp)),
                                 start(periranLocGR_200bp)-1, end(periranLocGR_200bp))
colnames(periranLoc_200bp_0based) <- c("chr", "start", "end")
write.table(periranLoc_200bp_0based, file = paste0(outDir, "periranLoc_200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)


