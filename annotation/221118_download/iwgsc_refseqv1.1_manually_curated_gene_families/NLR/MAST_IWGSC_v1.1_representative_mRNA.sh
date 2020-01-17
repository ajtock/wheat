#!/bin/bash

# Usage:
# ./MAST_IWGSC_v1.1_representative_mRNA.sh IWGSC_v1.1_nlr_representative_mRNA 

prefix=$1

mast /applications/NLR-Parser/NLR-Parser_1.0/meme.xml \
     ${prefix}_translate6frame.fasta \
     -o mast_out_${prefix} -ev 10000000
