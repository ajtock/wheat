#!/applications/R/R-3.5.0/bin/Rscript

# Compile alignment stats for multiple data sets into one supplemental table

ChIPNames <- c(
               "DMC1_Rep1_ChIP",
               "ASY1_CS_Rep1_ChIP",
               "H3K4me3_Rep1_ChIP",
               "H3K9me2_Rep1_ChIP",
               "H3K27me1_Rep1_ChIP",
               "H2AZ_Rep1_ChIP",
               "MNase_Rep1",
               "H3K9ac_ChIP_SRR6350667",
               "H3K27me3_ChIP_SRR6350666",
               "H3K36me3_ChIP_SRR6350670",
               "H3_input_SRR6350669",
               "H3K4me1_Rep1_ChIP_SRR8126618",
               "H3K27ac_Rep1_ChIP_SRR8126621",
               "CENH3_ChIP_SRR1686799"
              )
ChIPaccessions <- c(
                    "TBD",
                    "TBD",
                    "TBD",
                    "TBD",
                    "TBD",
                    "TBD",
                    "TBD",
                    "SRR6350667",
                    "SRR6350666",
                    "SRR6350670",
                    "SRR6350669",
                    "SRR8126618",
                    "SRR8126621",
                    "SRR1686799"
                   )
ChIPNamesDir <- c(
                  "DMC1",
                  "ASY1_CS",
                  "H3K4me3",
                  "H3K9me2",
                  "H3K27me1",
                  "H2AZ",
                  "MNase",
                  "H3K9ac",
                  "H3K27me3",
                  "H3K36me3",
                  "input",
                  "H3K4me1",
                  "H3K27ac",
                  "CENH3"
                 )
ChIPNamesTable <- c( 
                    "DMC1",
                    "ASY1",
                    "H3K4me3",
                    "H3K9me2",
                    "H3K27me1",
                    "H2A.Z",
                    "MNase",
                    "H3K9ac",
                    "H3K27me3",
                    "H3K36me3",
                    "Input",
                    "H3K4me1",
                    "H3K27ac",
                    "CENH3"
                   )
ChIPcite <- c(
              "This study",
              "This study",
              "This study",
              "This study",
              "This study",
              "This study",
              "This study",
              "Appels et. al. 2018 Science 361. doi: 10.1126/science.aar7191",
              "Appels et. al. 2018 Science 361. doi: 10.1126/science.aar7191",
              "Appels et. al. 2018 Science 361. doi: 10.1126/science.aar7191",
              "Appels et. al. 2018 Science 361. doi: 10.1126/science.aar7191",
              "Li et al. 2019 Genome Biol. 20. doi: 10.1186/s13059-019-1746-8",
              "Li et al. 2019 Genome Biol. 20. doi: 10.1186/s13059-019-1746-8",
              "Guo et al. 2016 PLoS Genet. 12. doi: 10.1371/journal.pgen.1005997"
             )
ChIPPEreadlen <- c( 
                   75,
                   75,
                   75,
                   75,
                   75,
                   75,
                   75,
                   75,
                   75,
                   75,
                   75,
                   150,
                   150,
                   100
                  )

ChIPStats <- sapply(seq_along(ChIPNames), function(x) {
  if(ChIPNames[x] %in% c("H3K4me3_ChIP_SRR6350668",
                         "H3K27me3_ChIP_SRR6350666",
                         "H3K36me3_ChIP_SRR6350670",
                         "H3K9ac_ChIP_SRR6350667",
                         "H3_input_SRR6350669",
                         "CENH3_ChIP_SRR1686799")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/logs/alignment_stats/", ChIPNames[x], "_alignment_summary.txt")
  } else if(ChIPNames[x] %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
                                "H3K27ac_Rep1_ChIP_SRR8126621")) {
    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/logs/alignment_stats/", ChIPNames[x], "_alignment_summary.txt")
  } else {
    paste0("/home/ajt200/analysis/wheat/",
           ChIPNamesDir[x], "/snakemake_ChIPseq/logs/alignment_stats/", ChIPNames[x], "_alignment_summary.txt")
  }
})

# Get genome size for calculating average depth of coverage
chrLens <- read.table("/home/ajt200/analysis/wheat/H3K4me3/snakemake_ChIPseq/data/index/wheat_v1.0.fa.sizes")[,2]
genomeSize <- sum(chrLens)

alnStatsDFlist <- lapply(seq_along(ChIPStats), function(x) {
  data.frame(Library = ChIPNamesTable[x],
             read.table(ChIPStats[x], header = T),
             Study = ChIPcite[x],
             Accession = ChIPaccessions[x],
             Read_length = ChIPPEreadlen[x],
             Read_lengthx2 = paste0("2x", ChIPPEreadlen[x], " bp"),
             stringsAsFactors = F)
})
alnStatsDF <- do.call(rbind, alnStatsDFlist)
alnStatsDF <- data.frame(alnStatsDF,
                         Alignment_rate = round(alnStatsDF$Primary_alignments/alnStatsDF$Trimmed,
                                                digits = 2) * 100,
                         Proportion_unique = round(alnStatsDF$MAPQ..23_mismatches..6 /
                                                   alnStatsDF$MAPQ..2_mismatches..6,
                                                   digits = 2) * 100,
                         Fold_coverage = round( ( alnStatsDF$MAPQ..2_mismatches..6 *
                                                  ( alnStatsDF$Read_length * 2 ) ) /
                                                    genomeSize,
                                               digits = 2),
                         stringsAsFactors = F)
colnames(alnStatsDF) <- c("Library", "Total sequenced read pairs",
                          "Deduplicated", "Trimmed",
                          "Primary alignments",
                          "Multi and unique alignments (MAPQ >= 2, mismatches <= 6)",
                          "Unique alignments (MAPQ >= 23, mismatches <= 6)",
                          "Study", "Accession",
                          "Read length", "Read length x2",
                          "Alignment rate (%)",
                          "Proportion unique (%)",
                          "Fold coverage")
write.table(alnStatsDF,
            file = "Supplemental_TableS1_alignment_summary.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(alnStatsDF,
          file = "Supplemental_TableS1_alignment_summary.csv",
          row.names = F, quote = F)
