#!/bin/bash

source activate ChIPseq_mapping
snakemake -p --cores 24
source deactivate
