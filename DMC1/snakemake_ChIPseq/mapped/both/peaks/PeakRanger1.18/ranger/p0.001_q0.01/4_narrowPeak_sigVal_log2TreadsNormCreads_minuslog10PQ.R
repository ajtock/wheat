#!/applications/R/R-3.3.2/bin/Rscript

# Adjust "signalValue" (region treatment reads; column 7) by region control reads
# signalValue = log2((region treatment reads+1)/(region control reads+1))
# Use this signalValue with caution - unknown whether these read counts are
# normalised by library size by PeakRanger
# May be preferrable to use qValue (FDR) or pValue in downstream analyses,
# including IDR calculation

# -log10 transform p-values and q-values

# Usage:
# ./4_narrowPeak_sigVal_log2TreadsNormCreads_minuslog10PQ.R DMC1_Rep1_ChIP 0.001 0.01

args <- commandArgs(trailingOnly = TRUE)
ChIPLibName <- args[1]
pval <- as.character(args[2])
qval <- as.character(args[3])

peaks <- read.table(paste0(ChIPLibName, "_peaks_peakranger_ranger_p",
                           pval, "_q", qval,
                           "_Treads_Creads.narrowPeak.UntransformedPQ"))
colnames(peaks) <- c("chr", "start0based", "end",
                     "name", "score", "strand",
                     "signalVal", "pValUntrans", "qValUntrans",
                     "summit0based", "treads", "creads")
peaks <- cbind(peaks[,1:6], log2((peaks[,11]+1)/(peaks[,12]+1)),
               -log10(peaks[,8]), -log10(peaks[,9]), peaks[,10])
colnames(peaks) <- c("chr", "start0based", "end",
                     "name", "score", "strand",
                     "TreadsNormCreads", "pVal", "qVal",
                     "summit0based")
peaks$pVal[which(!is.finite(peaks$pVal))] <- 323
peaks$qVal[which(!is.finite(peaks$qVal))] <- 323
head(peaks)
write.table(peaks, file = paste0(ChIPLibName, "_peaks_peakranger_ranger_p",
                                 pval, "_q", qval,
                                 "_log2TreadsNormCreads.narrowPeak"),
            col.names = F, row.names = F, sep = "\t", quote = F)
