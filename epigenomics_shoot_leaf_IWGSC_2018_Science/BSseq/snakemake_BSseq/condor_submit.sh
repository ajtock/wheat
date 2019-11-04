#!/bin/bash

source activate BSseq_mapping
snakemake -p --cores 32
conda deactivate
