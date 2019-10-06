#!/bin/bash

# Concatenate the Triticum aestivum (cv. Chinese Spring) nuclear (IWGSC refseq v1.0),
# chloroplast, and mitochondrial genomes, and the control sequences described in
# README_chloroplast_Lambda_pUC19_control_sequences.txt
# to create one bismark-indexable reference genome (wheat_v1.0_incl_organelles_controls.fa)
# NOTE: FASTA header lines have been shorted for chloroplast, mitochrondrion, lambda and pUC19
# in wheat_v1.0_incl_organelles_controls.fa

faDir=separate_nuclear_organelles_controls

cat $faDir/wheat_v1.0.fa \
    $faDir/wheat_CS_chloroplast_genome.fa \
    $faDir/wheat_CS_mitochondrial_genome.fa \
    $faDir/lambda_NEB.fa \
    $faDir/pUC19_NEB.fa \
    > wheat_v1.0_incl_organelles_controls.fa
