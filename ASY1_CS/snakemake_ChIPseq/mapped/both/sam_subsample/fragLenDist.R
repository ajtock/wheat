#!/applications/R/R-3.5.0/bin/Rscript

## Estimate fragment length size distribution using SAM column 9
## See https://www.biostars.org/p/65262/ for a discussion of this approach

# Usage
# ./fragLenDist.R ASY1_CS_Rep1_ChIP ASY1_CS_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_subsample.sam 21

args <- commandArgs(trailingOnly = T)
libName <- args[1]
samFile <- args[2]
ncol <- as.numeric(args[3])

#registerDoParallel(cores = length(libNames))
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())

# SAM format parsing in R - adapted from code in Michael Dondrup's post at https://www.biostars.org/p/932/
# adapt the "what" part to match your file, and adjust "ncol" to number of columns in SAM file
sam <- scan(file = samFile,
            what = rep(list(" "), ncol), fill = TRUE)
print("Mean fragment length:")
print(mean(abs(as.numeric(sam[[9]])), na.rm = T))

pdf(file = paste0(libName, "_fragLenDist_hist.pdf"), height = 5, width = 5)
par(mfrow = c(1, 1), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
hist(abs(as.numeric(sam[[9]])), breaks = 250, col = "grey60", border = NA, lwd = 2,
     xlab = paste(c(libName, " fragLen (bp)"), sep = ""), ylab = "Fragments", main = "", cex.lab = 2, cex.axis = 2)
abline(v = mean(abs(as.numeric(sam[[9]])), na.rm = T), col = "black", lty = 2, lwd = 1)
dev.off()

