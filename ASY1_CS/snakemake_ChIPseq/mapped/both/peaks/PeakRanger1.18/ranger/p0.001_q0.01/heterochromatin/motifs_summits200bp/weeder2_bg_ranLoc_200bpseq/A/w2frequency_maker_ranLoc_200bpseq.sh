#!/bin/bash

# Usage:
# ./w2frequency_maker_ranLoc_200bpseq.sh ASY1_CS_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_heterochromatin_summits200bp

prefix=$1
toolDir="/home/ajt200/tools/Weeder2.0"

[ -d FreqFiles ] || mkdir FreqFiles

$toolDir/w2frequency_maker ${prefix}_randomLoci.fa TA ds

mv *.freq FreqFiles
