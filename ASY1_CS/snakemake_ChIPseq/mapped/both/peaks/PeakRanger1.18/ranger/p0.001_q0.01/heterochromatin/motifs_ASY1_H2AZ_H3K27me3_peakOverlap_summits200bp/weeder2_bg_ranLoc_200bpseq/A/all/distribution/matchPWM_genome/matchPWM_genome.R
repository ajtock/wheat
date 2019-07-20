#!/applications/R/R-3.5.0/bin/Rscript

# Identify genomic loci that match motifs enriched at ASY1 peaks
# using motif position weight matrices (PWM)
# Profile base composition (proportion) in regions flanking these loci

# Usage:
# ./matchPWM_genome.R 'A' 'heterochromatin' 100 '100-bp'

#genomeName <- "A"
#region <- "heterochromatin"
#flankSize <- 100
#flankName <- "100-bp"

args <- commandArgs(trailingOnly = T)
genomeName <- args[1]
region <- args[2]
flankSize <- as.numeric(args[3])
flankName <- args[4]

library(Biostrings)
library(BSgenome.Taestivum.Cambridge.iwgscRefseqv1)
genome <- Taestivum
library(parallel)

plotDir <- "plots/"
propDir <- paste0(plotDir, "proportion/")
gffDir <- "gff/"
bedDir <- "bed/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
system(paste0("[ -d ", propDir, " ] || mkdir ", propDir))
system(paste0("[ -d ", gffDir, " ] || mkdir ", gffDir))
system(paste0("[ -d ", bedDir, " ] || mkdir ", bedDir))
plotDir <- propDir
 
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
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  regionGR <- regionGR[grep(genomeName,
                            seqnames(regionGR))@values]
} else if(region == "pericentromeres") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
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
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
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
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "pericentromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(genomeName,
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, heterochromatin, centromeres, pericentromeres, or genomewide")
}

# Extract genomeName-subgenome chromosomes
chrs <- chrs[grep(genomeName, chrs)]
genomeNums <- grep(genomeName,
                   names(genome))
chrList <- mclapply(genomeNums, function(i) {
  genome[[i]]
}, mc.cores = length(genomeNums))
names(chrList) <- names(genome)[genomeNums]

# Get number of motif PWMs
pwmNum <- as.numeric(system("ls -1 ../../MAT*.pwm | wc -l",
                            intern = T))

# Get PMM names
pwmName <- sapply(seq_along(1:pwmNum), function(x) {
  pwmName <- system(paste0("ls ../../MAT", x, "_*.pwm"),
                    intern = T) 
  pwmName <- substring(text = pwmName,
                       first = 7)
  sub(pattern = ".pwm",
      replacement = "",
      x = pwmName)
})

# Load PWMs as list
pwmList <- lapply(1:pwmNum, function(x) {
  as.matrix(read.table(system(paste0("ls ../../MAT", x, "_*.pwm"),
                              intern = T),
                       skip = 1, row.names = 1))
})

# Obtain motif matches in each genomeName-subgenome chromosome
mclapply(seq_along(pwmList), function(x) {
  hitsGR <- GRanges()
  for(i in seq_along(chrList)) {
    print(chrs[i])
    hits <- matchPWM(pwm = pwmList[[x]],
                     subject = chrList[[i]],
                     min.score = "99%")
    hitsGRchr <- GRanges(seqnames = chrs[i],
                         ranges = hits@ranges,
                         strand = "+")
    hitsRC <- matchPWM(pwm = reverseComplement(pwmList[[x]]),
                       subject = chrList[[i]],
                       min.score = "99%")
    hitsRCGRchr <- GRanges(seqnames = chrs[i],
                           ranges = hitsRC@ranges,
                           strand = "-")
    hitsGRchr <- sort(c(hitsGRchr, hitsRCGRchr))
    hitsGR <- append(hitsGR, hitsGRchr)
  }
  # Remove hits overlapping masked region
  mask_hits_overlap <- findOverlaps(query = maskGR,
                                    subject = hitsGR,
                                    type = "any",
                                    select = "all",
                                    ignore.strand = TRUE)
  hitsGR <- hitsGR[-subjectHits(mask_hits_overlap)]
  print(hitsGR)
  # Write as GRanges objects and in BED and GFF3 format
  save(hitsGR,
       file = paste0(pwmName[x], "_matchPWM_loci_in_",
                     genomeName, "genome_", region, ".RData"))
  hitsgff <- data.frame(chr = as.character(seqnames(hitsGR)),
                        source = as.character(rep(".")),
                        feature = as.character(rep(paste0(pwmName[x], "_match"))),
                        start = as.integer(start(hitsGR)),
                        end = as.integer(end(hitsGR)),
                        score = as.character(rep(".")),
                        strand = as.character(strand(hitsGR)),
                        phase = as.character(rep(".")),
                        attribute = as.character(rep(".")))
  write.table(hitsgff,
              file = paste0(gffDir, pwmName[x], "_matchPWM_loci_in_",
                            genomeName, "genome_", region, ".gff"),
              row.names = F, col.names = F, quote = F, sep = "\t")
  hitsbed <- data.frame(chr = as.character(seqnames(hitsGR)),
                        start = as.integer(start(hitsGR)-1),
                        end = as.integer(end(hitsGR)),
                        name = as.integer(1:length(hitsGR)),
                        score = rep("NA", length(hitsGR)),
                        strand = as.character(strand(hitsGR)))
  write.table(hitsbed,
              file = paste0(bedDir, pwmName[x], "_matchPWM_loci_in_",
                            genomeName, "genome_", region, ".bed"), 
              row.names = F, col.names = F, quote = F, sep = "\t")
}, mc.cores = length(pwmList))


# Load ASY1 peaks that overlap H2A.Z and H3K27me3 peaks
peaks <- read.table(paste0("/home/ajt200/analysis/wheat/ASY1_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
                           "ASY1_CS_Rep1_ChIP_rangerPeaksGR_minuslog10Qsorted_p0.001_q0.01_noMinWidth_in_",
                           genomeName, "genome_", region, ".gff"),
                    header = F)
peaksGR <- GRanges(seqnames = peaks$V1,
                   ranges = IRanges(start = peaks$V4,
                                    end = peaks$V5),
                   strand = "*")

# Load list of motif matches in GRanges format
hitsGRList <- lapply(seq_along(pwmName), function(x) {
  load(file = paste0(pwmName[x], "_matchPWM_loci_in_",
                     genomeName, "genome_", region, ".RData"))
})

library(doParallel)
# Change number of cores to reflect number of samples you want to process simultaneously
registerDoParallel(cores = length(pwmList))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(x = 1:pwmNum) %dopar% {
  print(x)
  load(file = paste0(pwmName[x], "_matchPWM_loci_in_",
                     genomeName, "genome_", region, ".RData"))
  peaks_hits_overlap <- findOverlaps(query = peaksGR,
                                     subject = hitsGR,
                                     type = "any",
                                     select = "all",
                                     ignore.strand = TRUE)
  hitsGR_peaksGR <- hitsGR[unique(subjectHits(peaks_hits_overlap))]
  print(hitsGR)
  print(hitsGR_peaksGR)
  hitsGRseq <- DNAStringSet()
  hitsGR_peaksGRseq <- DNAStringSet()
  for(h in seq_along(chrs)) {
    hitsGRseq_chr <- DNAStringSet()
    hitsGR_peaksGRseq_chr <- DNAStringSet()
    # Obtain sequences for each motif match and flanking flankSize-bp regions as DNAStringSet object
    # hits not overlapping peaks
    for(i in seq_along(hitsGR[seqnames(hitsGR) == chrs[h]])) {
      hitsGRseq_chr <- c(hitsGRseq_chr,
                         DNAStringSet(chrList[[h]][(start(ranges(hitsGR)[seqnames(hitsGR)
                                                          == chrs[h]][i])-flankSize):
                                                   (end(ranges(hitsGR)[seqnames(hitsGR)
                                                        == chrs[h]][i])+flankSize)]))
    }
    hitsGRseq <- c(hitsGRseq, hitsGRseq_chr)
    # hits overlapping peaks
    for(i in seq_along(hitsGR_peaksGR[seqnames(hitsGR_peaksGR) == chrs[h]])) {
      hitsGR_peaksGRseq_chr <- c(hitsGR_peaksGRseq_chr,
                                    DNAStringSet(chrList[[h]][(start(ranges(hitsGR_peaksGR)[seqnames(hitsGR_peaksGR)
                                                                     == chrs[h]][i])-flankSize):
                                                              (end(ranges(hitsGR_peaksGR)[seqnames(hitsGR_peaksGR)
                                                                   == chrs[h]][i])+flankSize)]))
    }
    hitsGR_peaksGRseq <- c(hitsGR_peaksGRseq, hitsGR_peaksGRseq_chr)
  }
  # Generate position frequency matrix (PFM)
  hitsGRpfm <- consensusMatrix(hitsGRseq)
  hitsGR_peaksGRpfm <- consensusMatrix(hitsGR_peaksGRseq)
  # Convert frequencies to proportions and retain rows 1:4
  hitsGRprm <- prop.table(hitsGRpfm, 2)[1:4,]
  hitsGR_peaksGRprm <- prop.table(hitsGR_peaksGRpfm, 2)[1:4,]
  # Re-order rows for stack barplot representation
  hitsGRprm_AGTC <- rbind(hitsGRprm[1,],
                          hitsGRprm[3,],
                          hitsGRprm[4,],
                          hitsGRprm[2,])
  rownames(hitsGRprm_AGTC) <- c("A", "G", "T", "C")
  hitsGR_peaksGRprm_AGTC <- rbind(hitsGR_peaksGRprm[1,],
                                  hitsGR_peaksGRprm[3,],
                                  hitsGR_peaksGRprm[4,],
                                  hitsGR_peaksGRprm[2,])
  rownames(hitsGR_peaksGRprm_AGTC) <- c("A", "G", "T", "C")
  pdf(paste0(plotDir,
             "ASY1_CS_Rep1_ChIP_peak_motif_", pwmName[x],
             #paste0(strsplit(consensusString(hitsGRseq),
             #                split = "")[[1]][(flankSize+1):(flankSize+1+mean(width(hitsGR))-1)],
             #       collapse = ""),
             "_base_proportions_in_", genomeName, "genome_", region, ".pdf"),
      height = 5, width = 20)
  par(mfrow = c(1, 2),
      mgp = c(1.5, 0.5, 0),
      mar = c(3.5, 3.1, 3.5, 2.1))
  barplot(hitsGR_peaksGRprm_AGTC,
          col = c("darkgreen", "orange", "red", "dodgerblue"),
          xlab = paste0("ASY1 peak motif ", pwmName[x],
                        #paste0(strsplit(consensusString(hitsGR_peaksGRseq),
                        #                split = "")[[1]][(flankSize+1):(flankSize+1+mean(width(hitsGR_peaksGR))-1)],
                        #       collapse = ""),
                        " matches and ", flankName, " flanks"),
          ylab = "Proportion",
          main = bquote("Matches overlapping ASY1 peaks in" ~
                        .(genomeName) * "-genome" ~ .(region) ~
                        "(" * italic("n") ~ "=" ~
                        .(prettyNum(length(hitsGR_peaksGR),
                                    big.mark = ",", trim = T)) *
                        ")"),
          legend.text = TRUE,
          args.legend = list(
            x = ncol(hitsGR_peaksGRprm_AGTC) + 60,
            y = 0.2,
            bty = "n"
          )
         )
  barplot(hitsGRprm_AGTC,
          col = c("darkgreen", "orange", "red", "dodgerblue"),
          xlab = paste0("ASY1 peak motif ", pwmName[x],
                        #paste0(strsplit(consensusString(hitsGR_peaksGRseq),
                        #                split = "")[[1]][(flankSize+1):(flankSize+1+mean(width(hitsGR_peaksGR))-1)],
                        #       collapse = ""),
                        " matches and ", flankName, " flanks"),
          ylab = "Proportion",
          main = bquote("All matches in" ~
                        .(genomeName) * "-genome" ~ .(region) ~
                        "(" * italic("n") ~ "=" ~
                        .(prettyNum(length(hitsGR),
                                    big.mark = ",", trim = T)) *
                        ")"),
          legend.text = TRUE,
          args.legend = list(
            x = ncol(hitsGRprm_AGTC) + 60,
            y = 0.2,
            bty = "n"
          )
         )
  dev.off()
}

