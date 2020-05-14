#!/applications/R/R-3.5.0/bin/Rscript

# Create a table containing wheat chromosome physical and genetic lengths (Mb and cM)

map <- read.table("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_mapping_data.txt",
                  header = T)
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]

chrLenTable <- data.frame()
for(x in 1:length(chrs)) {
  chrLenTable_chrx <- data.frame(Chromosome = chrs[x],
                                 Physical_length_in_Mb = as.numeric(chrLens[x]/1e6),
                                 Genetic_length_in_cM = max(map[map$chromosome == chrs[x],]$geneticPosition,
                                                            na.rm = T),
                                 # Or
                                 #Genetic_length_in_cM = map[map$chromosome == chrs[x],][
                                 #                        dim(map[map$chromosome == chrs[x],])[1],
                                 #                        dim(map[map$chromosome == chrs[x],])[2]],
                                 stringsAsFactors = F)
  chrLenTable <- rbind(chrLenTable, chrLenTable_chrx)
}
write.table(chrLenTable,
            file = "wheat_chromosome_physical_and_genetic_lengths.tsv",
            sep = "\t", quote = F, row.names = F)
