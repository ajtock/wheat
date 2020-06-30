#!/applications/R/R-3.5.0/bin/Rscript

# Extract 266 NLR-encoding genes that were upregulated in response to disease treatment
# (i.e., those with a log2-fold-increase >= 4)

library(rtracklayer)

# Load NLRs
NLRs <- read.table("NLR_genes_complete_representative_mRNA.gff3",
                   header = F, stringsAsFactors = F)
NLRs$V9 <- sub(pattern = "\\.\\d+", replacement = "",
               x = NLRs$V9)
NLRs$V9 <- sub(pattern = "2G", replacement = "1G",
               x = NLRs$V9)

# Load gene expression table
geneexp <- readRDS("/home/ajt200/analysis/wheat/RNAseq_RamirezGonzalez_Uauy_2018_Science/EI_grassroots_data_repo/TablesForExploration/MeanTpms.rds")
geneexp <- geneexp[!grepl("LC", geneexp$gene),]
geneexp <- geneexp[!grepl("CSU", geneexp$gene),]

# Get NLR-encoding genes upregulated upon challenge with
# powdery mildew
pm <- geneexp[geneexp$factor == "powdery mildew" &
              geneexp$subset == "disease_stress",]
pmc <- geneexp[geneexp$factor == "mildew stripe rust 2 control" &
               geneexp$subset == "disease_stress_control",]
stopifnot(identical(pm$gene, pmc$gene))
# Add offset where value is 0 to avoid infinite log2 ratios
# Order of these 4 commands is important:
#pm[which(pmc$value == 0),]$value <- pm[which(pmc$value == 0),]$value+1
#pmc[which(pmc$value == 0),]$value <- pmc[which(pmc$value == 0),]$value+1
#pmc[which(pm$value == 0),]$value <- pmc[which(pm$value == 0),]$value+1
#pm[which(pm$value == 0),]$value <- pm[which(pm$value == 0),]$value+1
pm_upreg <- pm[which(log2((pm$value)/(pmc$value)) >= 4),]
NLR_pm_upreg <- pm_upreg[which(pm_upreg$gene %in% NLRs$V9),]$gene

# stripe rust
sr <- geneexp[geneexp$factor == "stripe rust" &
              geneexp$subset == "disease_stress",]
src <- geneexp[geneexp$factor == "stripe rust control" &
               geneexp$subset == "disease_stress_control",]
stopifnot(identical(sr$gene, src$gene))
# Add offset where value is 0 to avoid infinite log2 ratios
# Order of these 4 commands is important:
#sr[which(src$value == 0),]$value <- sr[which(src$value == 0),]$value+1
#src[which(src$value == 0),]$value <- src[which(src$value == 0),]$value+1
#src[which(sr$value == 0),]$value <- src[which(sr$value == 0),]$value+1
#sr[which(sr$value == 0),]$value <- sr[which(sr$value == 0),]$value+1
sr_upreg <- sr[which(log2((sr$value)/(src$value)) >= 4),]
NLR_sr_upreg <- sr_upreg[which(sr_upreg$gene %in% NLRs$V9),]$gene

# stripe rust 2
sr2 <- geneexp[geneexp$factor == "stripe rust 2" &
               geneexp$subset == "disease_stress",]
sr2c <- geneexp[geneexp$factor == "mildew stripe rust 2 control" &
                geneexp$subset == "disease_stress_control",]
stopifnot(identical(sr2$gene, sr2c$gene))
# Add offset where value is 0 to avoid infinite log2 ratios
# Order of these 4 commands is important:
#sr2[which(sr2c$value == 0),]$value <- sr2[which(sr2c$value == 0),]$value+1
#sr2c[which(sr2c$value == 0),]$value <- sr2c[which(sr2c$value == 0),]$value+1
#sr2c[which(sr2$value == 0),]$value <- sr2c[which(sr2$value == 0),]$value+1
#sr2[which(sr2$value == 0),]$value <- sr2[which(sr2$value == 0),]$value+1
sr2_upreg <- sr2[which(log2((sr2$value)/(sr2c$value)) >= 4),]
NLR_sr2_upreg <- sr2_upreg[which(sr2_upreg$gene %in% NLRs$V9),]$gene

# Fusarium pseudograminearum
fp <- geneexp[geneexp$factor == "Fusarium pseudograminearum" &
              geneexp$subset == "disease_stress",]
fpc <- geneexp[geneexp$factor == "Fusarium pseudograminearum control" &
               geneexp$subset == "disease_stress_control",]
stopifnot(identical(fp$gene, fpc$gene))
# Add offset where value is 0 to avoid infinite log2 ratios
# Order of these 4 commands is important:
#fp[which(fpc$value == 0),]$value <- fp[which(fpc$value == 0),]$value+1
#fpc[which(fpc$value == 0),]$value <- fpc[which(fpc$value == 0),]$value+1
#fpc[which(fp$value == 0),]$value <- fpc[which(fp$value == 0),]$value+1
#fp[which(fp$value == 0),]$value <- fp[which(fp$value == 0),]$value+1
fp_upreg <- fp[which(log2((fp$value)/(fpc$value)) >= 4),]
NLR_fp_upreg <- fp_upreg[which(fp_upreg$gene %in% NLRs$V9),]$gene

# Fusarium pseudograminearum study2
fp2 <- geneexp[geneexp$factor == "Fusarium pseudograminearum study2" &
               geneexp$subset == "disease_stress",]
fp2c <- geneexp[geneexp$factor == "Fusarium pseudograminearum control study2" &
                geneexp$subset == "disease_stress_control",]
stopifnot(identical(fp2$gene, fp2c$gene))
# Add offset where value is 0 to avoid infinite log2 ratios
# Order of these 4 commands is important:
#fp2[which(fp2c$value == 0),]$value <- fp2[which(fp2c$value == 0),]$value+1
#fp2c[which(fp2c$value == 0),]$value <- fp2c[which(fp2c$value == 0),]$value+1
#fp2c[which(fp2$value == 0),]$value <- fp2c[which(fp2$value == 0),]$value+1
#fp2[which(fp2$value == 0),]$value <- fp2[which(fp2$value == 0),]$value+1
fp2_upreg <- fp2[which(log2((fp2$value)/(fp2c$value)) >= 4),]
NLR_fp2_upreg <- fp2_upreg[which(fp2_upreg$gene %in% NLRs$V9),]$gene

# Zymoseptoria tritici
zt <- geneexp[geneexp$factor == "Zymoseptoria tritici" &
              geneexp$subset == "disease_stress",]
ztc <- geneexp[geneexp$factor == "Zymoseptoria tritici control" &
               geneexp$subset == "disease_stress_control",]
stopifnot(identical(zt$gene, ztc$gene))
# Add offset where value is 0 to avoid infinite log2 ratios
# Order of these 4 commands is important:
#zt[which(ztc$value == 0),]$value <- zt[which(ztc$value == 0),]$value+1
#ztc[which(ztc$value == 0),]$value <- ztc[which(ztc$value == 0),]$value+1
#ztc[which(zt$value == 0),]$value <- ztc[which(zt$value == 0),]$value+1
#zt[which(zt$value == 0),]$value <- zt[which(zt$value == 0),]$value+1
zt_upreg <- zt[which(log2((zt$value)/(ztc$value)) >= 4),]
NLR_zt_upreg <- zt_upreg[which(zt_upreg$gene %in% NLRs$V9),]$gene

NLR_upreg_IDs <- sort(unique(c(NLR_pm_upreg, NLR_sr_upreg, NLR_sr2_upreg,
                               NLR_fp_upreg, NLR_fp2_upreg, NLR_zt_upreg)))
NLR_upreg_IDs <- sub("1G", "2G", NLR_upreg_IDs)

# Load representative genes
inDir <- "/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/"
mRNA_rep <- readGFF(paste0(inDir, "IWGSC_v1.1_HC_20170706_representative_mRNA.gff3"))
print(dim(mRNA_rep))
#[1] 107891      9
mRNA_rep$groupTrunc <- sub(pattern = "\\.\\d+", replacement = "",
                           x = mRNA_rep$group)

NLR_upreg <- mRNA_rep[mRNA_rep$groupTrunc %in% NLR_upreg_IDs,]
NLR_upreg <- NLR_upreg[order(NLR_upreg$seqid,
                             NLR_upreg$start,
                             NLR_upreg$end),]
print(dim(NLR_upreg))
write.table(NLR_upreg[,1:9],
            file = "NLR_genes_upregulated_disease_treatment_representative_mRNA.gff3",
            quote = F, sep = "\t", row.names = F, col.names = F)
NLR_upreg_bed <- data.frame(chr = as.character(NLR_upreg[,1]),
                            start = as.integer(NLR_upreg[,4]-1),
                            end = as.integer(NLR_upreg[,5]),
                            name = as.character(NLR_upreg[,9]),
                            score = as.numeric(NLR_upreg[,6]),
                            strand = as.character(NLR_upreg[,7]),
                            stringsAsFactors = F)
write.table(NLR_upreg_bed,
            file = "NLR_genes_upregulated_disease_treatment_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)
