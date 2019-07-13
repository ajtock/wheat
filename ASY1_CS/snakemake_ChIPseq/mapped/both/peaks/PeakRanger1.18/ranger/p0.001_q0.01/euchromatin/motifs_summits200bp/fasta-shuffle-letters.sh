#!/bin/bash

# Usage:
# ./fasta-shuffle-letters.sh ASY1_CS_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp 

prefix=$1
toolDir=/home/ajt200/meme/bin

$toolDir/fasta-shuffle-letters ${prefix}.fa \
                               ${prefix}_shuffled.fa
$toolDir/fasta-shuffle-letters ${prefix}_randomLoci.fa \
                               ${prefix}_randomLoci_shuffled.fa
