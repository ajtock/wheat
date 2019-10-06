#!/bin/bash

# Create FASTA index and get chromosome sizes

# Usage:
# ./samtools_faidx_chr_sizes.sh wheat_v1.0_incl_organelles_controls

i=$1

samtools faidx ${i}.fa
cut -f1,2 ${i}.fa.fai > ${i}.fa.sizes
