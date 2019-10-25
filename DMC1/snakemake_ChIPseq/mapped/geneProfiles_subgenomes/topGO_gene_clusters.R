#!/applications/R/R-3.5.0/bin/Rscript

########################################
# Analyse genes for GO term enrichment #
# relative to genes in given subgenome #
########################################

# This script is based on a post by Avril Coghlan:
# http://avrilomics.blogspot.co.uk/2015/07/using-topgo-to-test-for-go-term.html

# Usage:
# /applications/R/R-3.5.0/bin/Rscript ./topGO_gene_clusters.R BP 0.1 'log2_DMC1_Rep1_ChIP_control_in_promoters' 1 4 'genomewide' 'A'

# Doesn't work with mclapply or dopar

#ont <- "BP"
#sigLevel <- 0.1
#clusterBy <- "log2_DMC1_Rep1_ChIP_control_in_promoters"
#clusterFirst <- 1
#clusterLast <- 4
#region <- "genomewide"
#genomeName <- "A"

#options(echo=TRUE) # if you want to see commands in output file

args <- commandArgs(trailingOnly = TRUE)
ont <- args[1]
sigLevel <- args[2]
clusterBy <- args[3]
clusterFirst <- as.integer(args[4])
clusterLast <- as.integer(args[5])
region <- args[6]
genomeName <- args[7]

suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))

targets <- lapply(seq_along(clusterFirst:clusterLast), function(x) {
  paste0("clusters_by_", clusterBy, "/",
         "cluster", as.character(x),
         "_of_", as.character(clusterLast),
         "_by_", clusterBy, "_of_genes_in_",
         genomeName, "genome_", region, ".txt")
})

# Read in GO annotations for genes to define "gene universe"
geneID2GO <- readMappings(file = paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.0_FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned_in_",
                                        genomeName, "genome_", region,
                                        "_GO_IDs_no_chrUn.tsv"))
geneUniverse <- names(geneID2GO)

genesetGO <- function(target) {
  # Define list of genes of interest; file should contain a single column of gene identifiers
  genesOfInterest <- as.character(read.table(target)$V1)
  
  # Specify where genes of interest appear in the gene universe vector
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  # Build GOdata object in topGO
  capture.output(GOdata <- new("topGOdata", description = "Clustered genes", ontology = ont,
                               allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO),
                 file="/dev/null")
  
  # Access list of genes of interest
  #sg <- sigGenes(GOdata)
  #print(str(sg))
  #print(numSigGenes(GOdata))
  
  # Run Fisher's exact tests to determine GO term enrichment
  capture.output(resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultWeight <- runTest(GOdata, algorithm = "weight", statistic = "fisher"),
                 file="/dev/null")
  capture.output(resultTopGO <- runTest(GOdata, algorithm = "weight01", statistic = "fisher"),
                 file="/dev/null")
  
  # Count number of results where weight01 gives a P-value <= sigLevel
  mySummary <- summary(attributes(resultTopGO)$score <= as.numeric(sigLevel))
  if(length(mySummary) > 2) {
    numSignif <- as.integer(mySummary[[3]])
  } else {
    numSignif <- as.integer(0)
  }
  if(numSignif < 2) {
    numSignif <- as.integer(2)
  }

  # List significant results and write to file
  capture.output(enrichRes <- GenTable(GOdata,
                                       classicFisher = resultClassic,
                                       elimFisher = resultElim,
                                       weightFisher = resultWeight,
                                       topGOFisher = resultTopGO,
                                       orderBy = "topGOFisher",
                                       ranksOf = "weightFisher",
                                       topNodes = numSignif),
                 file="/dev/null")
  
  # WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
  basename <- basename(target)
  len <- nchar(basename)
  basename <- substr(basename, 1, len-4)
  
  out_name <- paste(basename, "GO", ont, "enrichment.tsv", sep="_")
  folder <- paste0(dirname(target), "/GO")
  system(paste0("[ -d ", folder, " ] || mkdir ", folder))
  folder2 <- paste0(folder, "/", basename, "_GO_", ont)
  system(paste0("[ -d ", folder2, " ] || mkdir ", folder2))
  
  capture.output(write.table(enrichRes, file = file.path(folder, out_name), sep = "\t",
                             row.names = FALSE, col.names = TRUE, quote = FALSE),
                 file="/dev/null")
  
  # Visualise the positions of the top 5 statistically significant GO terms in the GO hierarchy
  out_name2 <- paste(basename, "GO", ont, "enrichment", sep="_")
  printGraph(GOdata, resultTopGO, firstSigNodes = 5,
             fn.prefix = file.path(folder, out_name2), useInfo = "all", pdfSW = TRUE)
  
  # Extract gene IDs annotated with significantly enriched GO terms
  myTerms <- enrichRes$GO.ID
  myGenes <- genesInTerm(GOdata, myTerms)
  for(i in 1:length(myTerms)) {
    myTerm <- myTerms[i]
    myGenesForTerm <- myGenes[myTerm][[1]]
    myFactor <- myGenesForTerm %in% genesOfInterest
    myGenesForTermT <- myGenesForTerm[myFactor == TRUE]
    myGenesForTermT <- paste(myGenesForTermT, collapse = ",")
    myGenesForTermT <- paste(myTerm, myGenesForTermT, sep = "\t")
    out_name3 <- paste0(basename, "_GO_", ont, "_enrichment_", myTerm, ".txt")  
    write(myGenesForTermT, file = file.path(folder2, out_name3))
  }
}

# Apply genesetGO() function to each target (differentially expressed genes file)
lapply(seq_along(targets), function(x) {
  genesetGO(target = targets[[x]])
})
