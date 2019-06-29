#!/bin/bash

## PeakRanger v1.18
## Use peakranger ranger to call narrow peaks in ChIP-seq data, using the input to control for background

# Usage:
# csmit -m 100G -c 48 "./peaks_peakranger_ranger.sh ASY1_CS ASY1_CS_Rep1_ChIP input H3_input_SRR6350669 0.05 0.05 180 48"

ChIP=$1
ChIPLibName=$2
control=$3
controlLibName=$4
pval=$5
qval=$6
fragLen=$7
threads=$8

ChIPDir="/home/ajt200/analysis/wheat/"$ChIP"/snakemake_ChIPseq/mapped/both"
controlDir="/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/"$control"/snakemake_ChIPseq/mapped/both"
rangerDir="/home/ajt200/tools/PeakRanger-1.18/bin"

$rangerDir/peakranger ranger -d $ChIPDir/$ChIPLibName"_MappedOn_wheat_v1.0_lowXM_both_sort.bam" \
                             -c $controlDir/$controlLibName"_MappedOn_wheat_v1.0_lowXM_both_sort.bam" \
                             --format bam -o $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval \
                             -p $pval -q $qval -l $fragLen --pad -t $threads --verbose
