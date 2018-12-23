#!/bin/bash

# Build index for wheat reference genome (Wheat_IWGSC_WGA_v1.0_pseudomolecules)
# using bowtie version 1.2.2

/home/ajt200/anaconda3/bin/bowtie-build --threads 48 161010_Chinese_Spring_v1.0_pseudomolecules.fasta wheat_v1.0
/home/ajt200/anaconda3/bin/bowtie-build --threads 48 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta wheat_v1.0_parts
