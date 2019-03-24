#!/bin/bash

source activate ChIPseq_mapping
snakemake -p --cores 1
conda deactivate
