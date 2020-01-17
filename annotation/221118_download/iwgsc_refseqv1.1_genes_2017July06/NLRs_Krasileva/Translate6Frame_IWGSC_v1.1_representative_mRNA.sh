#!/bin/bash

# Usage:
# ./Translate6Frame_IWGSC_v1.1_representative_mRNA.sh NB_ARC_gene_IDs_IWGSC_v1_Ksenia_Krasileva_representative_mRNA

prefix=$1

Translate6Frame.jar --input ${prefix}.fasta \
                    --output ${prefix}_translate6frame.fasta
