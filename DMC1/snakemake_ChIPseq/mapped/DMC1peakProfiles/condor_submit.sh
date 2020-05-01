#!/bin/bash

source activate ChIPseq_mapping
snakemake -p --cores 32
conda deactivate
