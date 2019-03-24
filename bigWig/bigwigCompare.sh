#!/bin/bash

# Use bigwigCompare within the deepTools suite to obtain the
# ratio of the number of reads per bin for two bigWig files
# First install deepTools:
# conda install -c bioconda deeptools

# Usage (the third argument can be, e.g., ratio, log2, subtract):
# bash bigwigCompare.sh CENH3_ChIP_SRR1686799_MappedOn_wheat_v1.0_lowXM_both_sort_norm H3_input_SRR6350669_MappedOn_wheat_v1.0_lowXM_both_sort_norm ratio 1000000

prefix1=$1
prefix2=$2
operation=$3
binSize=$4

bigwigCompare --bigwig1 ${prefix1}.bw \
              --bigwig2 ${prefix2}.bw \
              --operation $operation \
              --binSize $binSize \
              --numberOfProcessors max \
              --outFileName ${prefix1}_${prefix2}_${operation}.bw \
              --outFileFormat bigwig           
