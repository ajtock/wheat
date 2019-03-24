#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./TE_family_frequency_chrProfiles.R 1Mb 1000000

winName <- "1Mb"
winSize <- 1000000

args <- commandArgs(trailingOnly = T)
winName <- args[1]
winSize <- as.numeric(args[2])

library(rtracklayer)
library(parallel)
library(GenomicRanges)
library(doParallel)

inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/"
outDirSuperfams <- "/home/ajt200/analysis/wheat/chromosomeProfiles/TEs/superfamilies/"
outDirSubfams <- "/home/ajt200/analysis/wheat/chromosomeProfiles/TEs/subfamilies/"
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
superfamName <- c("Gypsy_LTR",
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
superfamList <- mclapply(seq_along(superfamCode), function(x)
  {
    TEsGFF[grep(superfamCode[x], TEsGFF$compo),]
  }, mc.cores = 48)

# Create a list of TE superfamilies, wherein each list element
# is a vector of TE subfamilies
subfamList <- mclapply(seq_along(superfamList), function(x)
  {
    # Extract the first subfamily specified in the
    # composition ("compo") column
    unique(gsub(" .*$", "", superfamList[[x]]$compo))
  }, mc.cores = 48)

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
#chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
#chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,2])
centromereEnd <- as.vector(read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11.txt")[,3])

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

registerDoParallel(cores = length(superfamList))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(h = 1:length(superfamList)) %dopar% {
  superfamProfile <- NULL
  for(i in 1:length(chrs)) {
    # Count superfamily TEs within windows
    chrTEs <- superfamList[[h]][superfamList[[h]]$seqid == chrs[i],]
    chrTEsGR <- GRanges(seqnames = chrs[i],
                        ranges = IRanges(start = chrTEs$start,
                                         end = chrTEs$end),
                        strand = chrTEs$strand)
    chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
    winTEs <- countOverlaps(chrWindowsGR,
                            chrTEsGR,
                            ignore.strand = T)
    chrProfile <- data.frame(chr = as.character(chrs[i]),
                             window = as.integer(start(chrWindowsGR)),
                             winTEs = as.integer(winTEs))
    superfamProfile <- rbind(superfamProfile, chrProfile)
  }
  write.table(superfamProfile,
              file = paste0(outDirSuperfams,
                            "TE_frequency_per_", winName,
                            "_superfamily_", superfamName[h], "_", superfamCode[h],
                            ".txt"))

  for(j in 1:length(subfamList[[h]])) {
    print(h)
    print(j)
    print(subfamList[[h]][j])
    subfamProfile <- NULL
    for(i in 1:length(chrs)) {
    print(chrs[i])
    # Count subfamily TEs within windows
      chrTEs <- superfamList[[h]][superfamList[[h]]$seqid == chrs[i],]
      chrTEs_subfam <- chrTEs[grep(subfamList[[h]][j], chrTEs$compo),]
      chrWindowsGR <- windowsGR[seqnames(windowsGR) == chrs[i]]
      if(dim(chrTEs_subfam)[1] != 0) {
        chrTEs_subfamGR <- GRanges(seqnames = chrs[i],
                                   ranges = IRanges(start = chrTEs_subfam$start,
                                                    end = chrTEs_subfam$end),
                                   strand = chrTEs_subfam$strand)
        winTEs_subfam <- countOverlaps(chrWindowsGR,
                                       chrTEs_subfamGR,
                                       ignore.strand = T)
        chrProfile_subfam <- data.frame(chr = as.character(chrs[i]),
                                        window = as.integer(start(chrWindowsGR)),
                                        winTEs = as.integer(winTEs_subfam))
      } else {
        chrProfile_subfam <- data.frame(chr = as.character(chrs[i]),
                                        window = as.integer(start(chrWindowsGR)),
                                        winTEs = as.integer(0))
      }
      subfamProfile <- rbind(subfamProfile, chrProfile_subfam)
    }
    write.table(subfamProfile,
                file = paste0(outDirSubfams,
                              "TE_frequency_per_", winName,
                              "_superfamily_", superfamName[h], "_", superfamCode[h],
                              "_subfamily_", subfamList[[h]][j], 
                              ".txt"))
  }
}

