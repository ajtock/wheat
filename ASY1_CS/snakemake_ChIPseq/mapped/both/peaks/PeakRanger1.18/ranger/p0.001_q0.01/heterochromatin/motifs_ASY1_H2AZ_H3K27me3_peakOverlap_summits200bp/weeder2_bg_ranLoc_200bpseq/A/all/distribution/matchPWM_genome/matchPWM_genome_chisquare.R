# Load genomic loci that match motifs enriched at REC8 peaks and define random loci
# Profile base composition (-log10(chi-square P)) in regions flanking these loci

library(Biostrings)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(segmentSeq)
library(regioneR)
#library(zoo)
library(TTR)

inDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_HA_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/"
PFMDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_HA_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/PFM/"
plotDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_HA_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/plots/chisquare_P/"

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
seqlevels(genome) <- sub("Chr", "", seqlevels(genome))
mask <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))
seqlevels(mask) <- sub("Chr", "", seqlevels(mask))

chr1 <- Athaliana$Chr1
chr2 <- Athaliana$Chr2
chr3 <- Athaliana$Chr3
chr4 <- Athaliana$Chr4
chr5 <- Athaliana$Chr5
chr_list <- list()
chr_list[[1]] <- chr1
chr_list[[2]] <- chr2
chr_list[[3]] <- chr3
chr_list[[4]] <- chr4
chr_list[[5]] <- chr5

num_GR <- as.numeric(system(paste0("ls -l ", inDir, "motif*_matchPWM_GRanges.RData | wc -l"), intern = T))

library(doParallel)
# Change number of cores to reflect number of samples you want to process simultaneously
registerDoParallel(cores = num_GR)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(x = 1:num_GR) %dopar% {
  print(x)
  load(file = paste0(inDir, "motif", x, "_matchPWM_GRanges.RData"))
  ranLoc.GRanges <- randomizeRegions(motif.GRanges, genome = genome, per.chromosome = TRUE, allow.overlaps = TRUE)

  matchSeq <- DNAStringSet()
  ranLocSeq <- DNAStringSet()
  for(h in 1:5) {
    matchSeq.chr <- DNAStringSet()
    ranLocSeq.chr <- DNAStringSet()
    # Obtain sequences for each motif match, ranLoc and flanking 20-bp regions as DNAStringSet object
    #for(i in 1:10) {
    for(i in 1:length(motif.GRanges[seqnames(motif.GRanges) == h])) {
      matchSeq.chr <- c(matchSeq.chr, DNAStringSet(chr_list[[h]][(start(ranges(motif.GRanges)[seqnames(motif.GRanges) == h][i])-20):(end(ranges(motif.GRanges)[seqnames(motif.GRanges) == h][i])+20)]))
      ranLocSeq.chr <- c(ranLocSeq.chr, DNAStringSet(chr_list[[h]][(start(ranges(ranLoc.GRanges)[seqnames(ranLoc.GRanges) == h][i])-20):(end(ranges(ranLoc.GRanges)[seqnames(ranLoc.GRanges) == h][i])+20)]))
    }
    matchSeq <- c(matchSeq, matchSeq.chr)
    ranLocSeq <- c(ranLocSeq, ranLocSeq.chr)
  }
  # Generate position frequency matrix (PFM) for each sequence
  matchSeqPFM <- consensusMatrix(matchSeq)
  ranLocSeqPFM <- consensusMatrix(ranLocSeq)

  write.table(matchSeqPFM, file = paste0(PFMDir, "motif", x, "_matchSeqPFM_20bpFlank.txt"))
  write.table(ranLocSeqPFM, file = paste0(PFMDir, "motif", x, "_ranLocSeqPFM_20bpFlank.txt"))

  #X2 = Σ (observed – expected)2 / expected

  matchMaxBases <- NULL
  matchMaxFreq <- NULL
  matchSumFreq <- NULL
  ranLocBases <- NULL
  ranLocFreq <- NULL
  ranLocSumFreq <- NULL
  for(i in 1:dim(matchSeqPFM)[2]) {
    matchMaxBases <- c(matchMaxBases, names(matchSeqPFM[,i][matchSeqPFM[,i] == max(matchSeqPFM[,i])]))
    ranLocBases <- c(ranLocBases, names(ranLocSeqPFM[,i][matchSeqPFM[,i] == max(matchSeqPFM[,i])]))
    matchMaxFreq <- c(matchMaxFreq, as.integer(matchSeqPFM[,i][matchSeqPFM[,i] == max(matchSeqPFM[,i])]))
    ranLocFreq <- c(ranLocFreq, as.integer(ranLocSeqPFM[,i][matchSeqPFM[,i] == max(matchSeqPFM[,i])]))
    matchSumFreq <- c(matchSumFreq, sum(matchSeqPFM[,i]))
    ranLocSumFreq <- c(ranLocSumFreq, sum(ranLocSeqPFM[,i]))
  }
 
  observed <- NULL
  expectedProportions <- NULL
  for(j in 1:length(matchMaxFreq)) {
    observed <- rbind(observed, c(matchMaxFreq[j], matchSumFreq[j]-matchMaxFreq[j]))
    expectedProportions <- rbind(expectedProportions,
                                 c(ranLocFreq[j]/ranLocSumFreq[j],
                                  (ranLocSumFreq[j]-ranLocFreq[j])/ranLocSumFreq[j]))
  }

  chisqList <- lapply(seq_along(1:dim(observed)[1]), function(k) {
    chisq.test(observed[k,1:2], p = expectedProportions[k,1:2])
  })
  save(chisqList, file = paste0(plotDir, "motif", x, "_chisqList_20bpFlank.RData"))

  chisqMatrix <- NULL
  for(k in 1:dim(observed)[1]) {
    chisq <- chisq.test(observed[k,1:2], p = expectedProportions[k,1:2])
    chisqMatrix <- rbind(chisq, chisqMatrix)
  }
  chisqDataFrame <- as.data.frame(chisqMatrix, row.names = as.character(1:18))
  save(chisqDataFrame, file = paste0(plotDir, "motif", x, "_chisqDataFrame_20bpFlank.RData"))
}

