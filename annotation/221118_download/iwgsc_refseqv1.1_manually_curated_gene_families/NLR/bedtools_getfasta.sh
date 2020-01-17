#!/bin/bash

# Extract nucleotide sequence at given coordinates (specified in a BED file)
# in the genome in fasta format

# Note that BED files use 0-based half-open coordinates;
# start coordinates are 0-based and end coordinates are 1-based,
# such that a feature START coordinate that cooresponds to the first base in
# a chromosome is numbered 0, while a feature END coordinate that corresponds
# to the second base in a chromosome is numbered 2 in the BED file:
# see https://genome.ucsc.edu/FAQ/FAQformat.html#format1

# Usage:
# ./bedtools_getfasta.sh 161010_Chinese_Spring_v1.0_pseudomolecules.fasta IWGSC_v1.1_nlr_representative_mRNA 

genome=$1
prefix=$2

/home/ajt200/anaconda3/bin/bedtools getfasta -fi ${genome} \
                                             -bed ${prefix}.bed \
                                             -fo ${prefix}.fasta \
                                             -name+
