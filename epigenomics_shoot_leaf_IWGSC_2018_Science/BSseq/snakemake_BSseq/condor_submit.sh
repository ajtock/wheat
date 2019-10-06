#!/bin/bash

source activate BSseq_mapping
snakemake -p --cores 48
conda deactivate
