#!/applications/R/R-3.4.0/bin/Rscript

# Create GFF3 file of features sorted by coverage levels in promoters,
# within feature bodies, or in terminators

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_sorted_by_log2_ChIP_control_writeGFF3.R H2AZ H2AZ_Rep1_ChIP ASY1 ASY1_CS_Rep1_ChIP p0.001_q0.01 2000 2kb '2 kb' 20 20bp bodies 'A' 'euchromatin'

#markName <- "H2AZ"
#libName <- "H2AZ_Rep1_ChIP"
#featureName <- "ASY1"
#featureLibName <- "ASY1_CS_Rep1_ChIP"
#sigLevel <- "p0.001_q0.01"
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
markName <- args[1]
libName <- args[2]
featureName <- args[3]
featureLibName <- args[4]
sigLevel <- args[5]
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
flankNamePlot <- args[8]
binSize <- as.numeric(args[9])
binName <- args[10]
region <- args[11]
genomeName <- args[12]
chrRegion <- args[13]

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
           "input/snakemake_ChIPseq/mapped/", featureName, "peakProfiles/matrices/")
  } else if(libNameControl[x] == "MNase_Rep1") {
    paste0("/home/ajt200/analysis/wheat/",
           "MNase/snakemake_ChIPseq/mapped/", featureName, "peakProfiles/matrices/")
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
                              featureName, "_CS_peaks_in_", genomeName, "genome_", chrRegion,
                              "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(libNameControl))

## Load libName coverage matrix, log2-transform, and extract region for sorting of features
if(libName %in% c("H3K4me3_ChIP_SRR6350668",
                  "H3K27me3_ChIP_SRR6350666",
                  "H3K36me3_ChIP_SRR6350670",
                  "H3K9ac_ChIP_SRR6350667",
                  "CENH3_ChIP_SRR1686799")) {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
                       markName,
                       "/snakemake_ChIPseq/mapped/", featureName, "peakProfiles/matrices/")
} else {
  covDirChIP <- paste0("/home/ajt200/analysis/wheat/",
                       markName,
                       "/snakemake_ChIPseq/mapped/", featureName, "peakProfiles/matrices/")
}

mat1 <- as.matrix(read.table(paste0(covDirChIP,
                                    libName,
                                    "_MappedOn_wheat_v1.0_lowXM_both_sort_norm_",
                                    featureName, "_CS_peaks_in_", genomeName, "genome_", chrRegion,
                                    "_matrix_bin", binName,
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
peaksbed <- read.table(paste0("/home/ajt200/analysis/wheat/", featureName,
                              "_CS/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/", sigLevel, "/",
                              featureLibName, "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel, "_noMinWidth_in_",
                              genomeName, "genome_", chrRegion, ".bed"),
                       header = F)
colnames(peaksbed) <- c("chr", "start", "end", "name", "score", "strand")
peaksgff <- data.frame(chr = as.character(peaksbed$chr),
                       source = as.character(rep(".")),
                       feature = as.character(rep(paste0(featureLibName, "_peak"))),
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
            file = paste0(gffDir,
                          featureLibName, "_rangerPeaksGRmergedOverlaps_minuslog10_", sigLevel, "_noMinWidth_in_",
                          genomeName, "genome_", chrRegion,
                          "_ordered_by_log2_", libName, "_control_in_bodies.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")
