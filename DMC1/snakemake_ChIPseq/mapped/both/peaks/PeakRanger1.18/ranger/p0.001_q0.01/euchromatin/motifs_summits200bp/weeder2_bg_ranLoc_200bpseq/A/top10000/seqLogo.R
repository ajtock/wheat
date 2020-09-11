#!/applications/R/R-3.5.0/bin/Rscript

# Convert weeder2 position weight matrices (PWMs)
# (that have been split into separate files using split_pwm.py)
# into sequence logos using seqLogo Bioconductor package

# Usage:
# ./seqLogo.R DMC1_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp

args <- commandArgs(trailingOnly = T)
featurePrefix <- args[1]

library(seqLogo)

seqLogoDir <- "seqLogo/"
system(paste0("[ -d ", seqLogoDir, " ] || mkdir ", seqLogoDir))

# Get number of PWMs in current working directory
num_pwm <- as.numeric(system("ls -1 MAT*.pwm | wc -l",
                             intern = T))

# Create individual motif logo files in EPS format
for(i in 1:num_pwm) {
  pwm_name <- system(paste0("ls MAT", i, "_*.pwm"),
                     intern = T)
  pwm_name <- sub(pattern = ".pwm",
                  replacement = "",
                  x = pwm_name) 
  pwm <- read.table(system(paste0("ls MAT", i, "_*.pwm"),
                           intern = T),
                    skip = 1, row.names = 1)
  print(head(pwm))
  postscript(paste0(seqLogoDir,
                    "MAT", i, ".eps"))
  seqLogo(pwm, xfontsize = 60, yfontsize = 60)
  dev.off()
}

pwm_list <- lapply(1:num_pwm, function(i) {
  read.table(system(paste0("ls MAT", i, "_*.pwm"),
                    intern = T),
             skip = 1, row.names = 1)
})

# Edit seqLogo command to allow logos to be combined in one PDF
# See https://support.bioconductor.org/p/35240/
mySeqLogo <- seqLogo::seqLogo

bad <- (sapply( body(mySeqLogo), "==", "grid.newpage()") |
         sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
body(mySeqLogo)[bad] <- NULL

pdf(paste0(seqLogoDir, featurePrefix,
           "_weeder2_motif_seqLogos.pdf"), height = 5*length(pwm_list), width = 10)
# par() doesn't work here 
# par(mfrow = c(10, 1)),
#     mar = c(3, 3, 2, 2),
#     mgp = c(1, 0.5, 0))
grid.newpage()
for(i in 1:length(pwm_list)) {
  pushViewport(viewport(x = 0.5,
                        y = 1-(i/length(pwm_list))+((1/length(pwm_list))/2),
                        width = 1, height = 1/length(pwm_list)))
  mySeqLogo(pwm_list[[i]],
            xaxis = F, yaxis = T,
            xfontsize = 60, yfontsize = 60)
  grid.text(sprintf("Motif %d", i),
            x = 0.6, y = 0.9,
            gp = gpar(fontsize = 60))
  popViewport()
}
dev.off()

