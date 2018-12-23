#!/bin/bash

# Usage:
# cat_fastq_gz_gzip_best.sh 1575_LIB18613_LDI16228_GGCTAC R1 CS+_2_LIB18613_LDI16228_R1

inName=$1
readN=$2
outName=$3

if [ ! -f "$outName.fastq.gz" ]; then
  cat $inName"_L001_"$readN".fastq" \
      $inName"_L002_"$readN".fastq" \
      | gzip -c -k --best > $outName.fastq.gz;
else
  echo "skipping $outName"
fi

