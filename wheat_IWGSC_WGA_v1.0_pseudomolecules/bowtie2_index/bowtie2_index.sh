#!/bin/bash

# Build index for wheat reference genome (Wheat_IWGSC_WGA_v1.0_pseudomolecules)
# using bowtie2 version 2.3.4.3

genome=$1
idxBaseName=$2

/home/ajt200/anaconda3/bin/bowtie2-build --threads 48 $genome $idxBaseName
