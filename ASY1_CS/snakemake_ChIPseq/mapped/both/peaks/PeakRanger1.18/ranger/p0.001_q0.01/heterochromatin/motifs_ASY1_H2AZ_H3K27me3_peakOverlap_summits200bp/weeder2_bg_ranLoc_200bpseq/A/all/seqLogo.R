#!/applications/R/R-3.5.0/bin/Rscript

# Convert weeder2 position weight matrices (PWMs)
# (that have been split into separate files using split_pwm.py)
# into sequence logos using seqLogo Bioconductor package

# Usage:
# ./seqLogo.R ASY1_CS_Rep1_ChIP_peaks_ASY1_CS_Rep1_ChIPsorted_in_Agenome_heterochromatin_summits200bp

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
  pwm_name <- system(paste0(" ls MAT", i, "_*.pwm"),
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
  seqLogo(pwm)
  dev.off()
}

pwm_list <- lapply(1:num_pwm, function(i) {
  read.table(system(paste0("ls MAT", i, "_*.pwm"),
                    intern = T),
             skip = 1, row.names = 1)
})
pdf(paste0(seqLogoDir, featurePrefix,
           "_weeder2_motif_seqLogos.pdf"))
for(i in 1:num_pwm) {
  seqLogo(pwm_list[[i]])
}
dev.off()


