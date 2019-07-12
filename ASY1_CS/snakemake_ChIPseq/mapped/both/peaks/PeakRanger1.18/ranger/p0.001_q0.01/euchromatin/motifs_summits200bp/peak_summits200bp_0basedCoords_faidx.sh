#!/bin/bash

execPATH="/applications/anaconda/anaconda3/bin"

# obtain locus sequences in fasta format with coordinates in 0-based bed file

$execPATH/faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed armrangerPeaks_summits200bp_0based.bed --out armrangerPeaks_summits200bpseq.fa
$execPATH/faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed perirangerPeaks_summits200bp_0based.bed --out perirangerPeaks_summits200bpseq.fa

$execPATH/faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed armranLoc_200bp_0based.bed --out armranLoc_200bpseq.fa
$execPATH/faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed periranLoc_200bp_0based.bed --out periranLoc_200bpseq.fa

