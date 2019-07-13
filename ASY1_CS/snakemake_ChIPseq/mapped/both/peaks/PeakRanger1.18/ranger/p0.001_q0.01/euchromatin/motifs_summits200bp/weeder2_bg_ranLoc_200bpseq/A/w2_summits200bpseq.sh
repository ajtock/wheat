#!/bin/bash

# Usage:
# ./w2_summits200bpseq.sh ASY1_CS_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp
# ./w2_summits200bpseq.sh ASY1_CS_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp_shuffled

prefix=$1
toolDir="/home/ajt200/tools/Weeder2.0"

weeder2 -f ${prefix}.fa \
        -O TA -chipseq -top 10000
