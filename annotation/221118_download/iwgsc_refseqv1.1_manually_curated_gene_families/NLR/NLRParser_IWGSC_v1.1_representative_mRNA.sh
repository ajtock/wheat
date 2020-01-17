#!/bin/bash

# Usage:
# ./NLRParser_IWGSC_v1.1_representative_mRNA.sh IWGSC_v1.1_nlr_representative_mRNA

prefix=$1

NLR-Parser.jar --input mast_out_${prefix}/mast.xml \
               --output ${prefix}_NLRParser_out.tsv \
               --sequenceFile ${prefix}_translate6frame.fasta \
               --pValue 1E-5 \
               --splitpattern _frame
