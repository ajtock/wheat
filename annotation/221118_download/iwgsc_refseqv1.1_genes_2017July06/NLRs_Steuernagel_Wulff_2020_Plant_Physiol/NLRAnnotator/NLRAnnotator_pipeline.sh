#!/bin/bash

# Apply NLR-Annotator to identify NLR classes and domains within a given set of sequences
# Further details and download avaialble at:
# https://github.com/steuernb/NLR-Annotator

# Usage:
# ./NLRAnnotator_pipeline.sh NLR_genes_representative_mRNA

prefix=$1

ChopSequence.jar -i ${prefix}.fasta \
                 -o ${prefix}_subseqs.fasta 

NLR-Parser_2.3.jar -t 48 \
                   -y /home_old/srna/meme/bin/mast \
                   -x /applications/NLR-Annotator/NLR-Annotator_0.7/meme.xml \
                   -i ${prefix}_subseqs.fasta \
                   -c ${prefix}_subseqs_mast.xml

NLR-Annotator.jar -i ${prefix}_subseqs_mast.xml \
                  -o ${prefix}_NLRAnnotator.tsv \
                  -g ${prefix}_NLRAnnotator.gff \
                  -b ${prefix}_NLRAnnotator.bed \
                  -m ${prefix}_NLRAnnotator_motifs.bed \
                  -f ${prefix}.fasta ${prefix}_NLRAnnotator.fasta 0
