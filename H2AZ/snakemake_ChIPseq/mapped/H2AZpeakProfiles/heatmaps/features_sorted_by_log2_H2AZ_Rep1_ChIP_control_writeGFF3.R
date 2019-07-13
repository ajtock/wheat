#!/applications/R/R-3.4.0/bin/Rscript

# Create GFF3 file of features sorted by coverage levels in promoters,
# within feature bodies, or in terminators

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_sorted_by_log2_H2AZ_Rep1_ChIP_control_writeGFF3.R H2AZ_Rep1_ChIP H2AZ_peaks_in_Agenome_euchromatin 2000 2kb '2 kb' 20 20bp bodies 'A' 'euchromatin'

#libName <- "H2AZ_Rep1_ChIP"
#featureName <- "H2AZ_peaks_in_Agenome_euchromatin"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "bodies"
#genomeName <- "A"
#chrRegion <- "euchromatin"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
featureName <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
region <- args[8]
genomeName <- args[9]
chrRegion <- args[10]

library(parallel)

gffDir <- "gff/"
regiongffDir <- paste0(gffDir, region, "_by_log2_", libName, "_control/")
system(paste0("[ -d ", gffDir, " ] || mkdir ", gffDir))
system(paste0("[ -d ", regiongffDir, " ] || mkdir ", regiongffDir))
gffDir <- regiongffDir

## Load Control coverage matrices to be used for log2 transformation
libNameControl <- c(
                    "MNase_Rep1",
                    "H3_input_SRR6350669"
                   )
markControl <- c(
                 "MNase",
                 "input"
                )
covDirControl <- sapply(seq_along(libNameControl), function(x) {
  if(libNameControl[x] == "H3_input_SRR6350669") {
    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
           "input/snakemake_ChIPseq/mapped/H2AZpeakProfiles/matrices/")
  } else if(libNameControl[x] == "MNase_Rep1") {
    paste0("/home/ajt200/analysis/wheat/",
           "MNase/snakemake_ChIPseq/mapped/H2AZpeakProfiles/matrices/")
  } else {
    if(!(libNameControl %in% c("H3_input_SRR6350669", "MNase_Rep1"))) {
      stop("libNameControl[x] is neither H3_input_SRR6350669 nor MNase_Rep1")
    }
  }
})

covMatControl <- mclapply(seq_along(libNameControl), function(x) {
  as.matrix(read.table(paste0(covDirControl[x],
                              libNameControl[x],
                              "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                              featureName, "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(libNameControl))

## Load libName coverage matrix, log2-transform, and extract region for sorting of features
mat1 <- as.matrix(read.table(paste0("../matrices/",
                                    libName,
                                    "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                                    featureName, "_matrix_bin", binName,
                                    "_flank", flankName, ".tab"),
                             header = F, skip = 3))

# Calculate log2(ChIP/control) for mat1, to be used for sorting features
# by decreasing levels
if(libName %in% c("ASY1_CS_Rep1_ChIP",
                  "MNase_Rep1",
                  "H3K4me3_ChIP_SRR6350668",
                  "H3K27me3_ChIP_SRR6350666",
                  "H3K36me3_ChIP_SRR6350670",
                  "H3K9ac_ChIP_SRR6350667",
                  "CENH3_ChIP_SRR1686799")) {
  # log2(ChIP/input)
  log2mat1 <- log2((mat1+1)/(covMatControl[[2]]+1))
} else {
  # log2(ChIP/MNase)
  log2mat1 <- log2((mat1+1)/(covMatControl[[1]]+1))
}

bodyLength <- (dim(log2mat1)[2]-((upstream+downstream)/binSize))*binSize
if( region == "promoters" ) {
  log2mat1Region <- log2mat1[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  log2mat1Region <- log2mat1[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
} else if ( region == "bodies" ) {
  log2mat1Region <- log2mat1[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else {
  print("The region name provided does not match 'promoters',  'terminators', or 'bodies'")
}
log2mat1RegionRowMeans <- rowMeans(log2mat1Region, na.rm = T)
log2mat1RegionRowMeansSorted <- sort.int(log2mat1RegionRowMeans,
                                         decreasing = T,
                                         index.return = T,
                                         na.last = T)

# Load peaks in BED format, convert into GFF3 format,
# and sort by decreasing log2mat1RegionRowMeans
# Note addition of 1 to 0-based BED start coordinates
peaksbed <- read.table(paste0("/home/ajt200/analysis/wheat/H2AZ/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.05_q0.05/",
                              libName,
                              "_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth_in_",
                              genomeName, "genome_", chrRegion, ".bed"),
                       header = F)
colnames(peaksbed) <- c("chr", "start", "end", "name", "score", "strand")
peaksgff <- data.frame(chr = as.character(peaksbed$chr),
                       source = as.character(rep(".")),
                       feature = as.character(rep(paste0(libName, "_peak"))),
                       start = as.integer(peaksbed$start+1),
                       end = as.integer(peaksbed$end),
                       score = as.numeric(log2mat1RegionRowMeans),
                       strand = as.character(rep(".")),
                       frame = as.character(rep(".")),
                       attribute = as.character(rep(".")))
peaksgffSorted <- peaksgff[order(peaksgff$score,
                                 decreasing = T,
                                 na.last = T),]
write.table(peaksgffSorted,
            file = paste0(gffDir, libName,
                          "_rangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth_in_",
                          genomeName, "genome_", chrRegion,
                          "_ordered_by_log2_", libName, "_control_in_bodies.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")
