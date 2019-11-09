#!/bin/bash

# Usage:
# csmit -m 10G -c 1 "bash ./BSseq_zcat_fastq_gz_best.sh BSseq_PRJNA436361 R1"

prefix=$1
readNo=$2
outName=$prefix"_"$readNo

if [ ! -f "$outName.fastq.gz" ]; then
  zcat SRR6792673_$readNo".fastq.gz" \
       SRR6792674_$readNo".fastq.gz" \
       SRR6792675_$readNo".fastq.gz" \
       SRR6792676_$readNo".fastq.gz" \
       SRR6792677_$readNo".fastq.gz" \
       SRR6792678_$readNo".fastq.gz" \
       SRR6792679_$readNo".fastq.gz" \
       SRR6792680_$readNo".fastq.gz" \
       SRR6792681_$readNo".fastq.gz" \
       SRR6792682_$readNo".fastq.gz" \
       SRR6792683_$readNo".fastq.gz" \
       SRR6792684_$readNo".fastq.gz" \
       SRR6792685_$readNo".fastq.gz" \
       SRR6792686_$readNo".fastq.gz" \
       SRR6792687_$readNo".fastq.gz" \
       SRR6792688_$readNo".fastq.gz" \
       SRR6792689_$readNo".fastq.gz" \
       | gzip --stdout --keep --best > $outName.fastq.gz;
else
  echo "skipping $outName"
fi
