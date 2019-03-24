#!/bin/bash

# Convert PeakRanger peak loci details file to narrowPeak format

# Usage:
# ./detailsTOnarrowPeak_Treads_Creads.sh H3K4me3_Rep1_ChIP 0.05 0.05

ChIPLibName=$1
pval=$2
qval=$3

tail -n +24 $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_details" \
  > $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead"
grep -v 'chrUn' $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead" \
  > $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead_tmp1"
awk 'BEGIN {OFS="\t"}; {$1 = $1; print}' $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead_tmp1" \
  > $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead_tmp2"
awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $8, $8, $8, $8, $6, $7, $5, $9, $10}' $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead_tmp2" \
  > $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead_tmp3"
awk 'BEGIN {OFS="\t"}; {$4 = "."; $5 = "."; $6 = "."; $7 = "."; $10 = $10-$2; print}' $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead_tmp3" \
  > $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_Treads_Creads.narrowPeak.UntransformedPQ"
rm $ChIPLibName"_peaks_peakranger_ranger_p"$pval"_q"$qval"_noHead"*

