#!/bin/bash

# Apply NLR-Parser to identify NLR classes and domains within a given set of sequences
# Further details and download avaialble at:
# https://github.com/steuernb/NLR-Parser

# Usage:
# ./NLRParser_pipeline.sh NB_ARC_genes_IWGSC_v1_Ksenia_Krasileva_representative_mRNA

prefix=$1

Translate6Frame.jar --input ${prefix}.fasta \
                    --output ${prefix}_translate6frame.fasta

mast /applications/NLR-Parser/NLR-Parser_1.0/meme.xml \
     ${prefix}_translate6frame.fasta \
     -o mast_out_${prefix} -ev 1000000000000000000

NLR-Parser.jar --input mast_out_${prefix}/mast.xml \
               --output ${prefix}_NLRParser_out.tsv \
               --sequenceFile ${prefix}_translate6frame.fasta \
               --pValue 1E-5 \
               --splitpattern _frame
