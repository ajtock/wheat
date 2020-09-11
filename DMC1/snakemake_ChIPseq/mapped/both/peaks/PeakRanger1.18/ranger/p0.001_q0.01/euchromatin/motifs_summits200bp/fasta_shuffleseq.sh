#!/bin/bash

# Use EMBOSS:6.6.0.0 shuffleseq to shuffle fasta sequences, maintaining base composition

# Usage:
# ./fasta_shuffleseq.sh DMC1_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp 

prefix=$1
toolDir=/home/ajt200/tools/shuffleseq

$toolDir/shuffleseq -sequence ${prefix}.fa -outseq ${prefix}_shuffled.fa
#$toolDir/fasta-shuffle-letters ${prefix}_randomLoci.fa ${prefix}_randomLoci_shuffled.fa
