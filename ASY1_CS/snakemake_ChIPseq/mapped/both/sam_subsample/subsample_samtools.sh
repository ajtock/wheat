#!/bin/bash

# Subsample bam files

# samtools 1.9

# Usage:
# csmit -m 1G -c 1 "./subsample_samtools.sh ASY1_CS_Rep1_ChIP"

prefix=$1
echo $prefix

samtoolsDir=/home/ajt200/anaconda3/envs/ChIPseq_mapping/bin
bamDir=..

proportion=0.1
randomSeed=85363
seedPlusProportion=$(echo $randomSeed + $proportion | bc -l)
echo $seedPlusProportion

# Subsample
$samtoolsDir/samtools view -s $seedPlusProportion \
                           $bamDir/${prefix}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
  > ${prefix}_MappedOn_wheat_v1.0_lowXM_both_sort_subsample.sam
