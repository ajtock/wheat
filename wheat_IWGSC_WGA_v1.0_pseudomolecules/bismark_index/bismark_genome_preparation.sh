#!/bin/bash

# NOTE: This requires a lot of RAM (> 260GB) - must run on node10 (500GB)

# Use Bismark (v0.20.0) to run bisulfite conversion of reference genome
# (e.g., wheat_v1.0_incl_organelles_controls.fa)
# and create bisulfite reference indexes for bowtie2

# Usage:
# csmit -m 300G -c 2 "bash bismark_genome_preparation.sh" 

toolDir=/home/ajt200/anaconda3/envs/BSseq_mapping/bin

$toolDir/bismark_genome_preparation --path_to_bowtie $toolDir/ --bowtie2 \
                                    --genomic_composition \
                                    /home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/bismark_index/
