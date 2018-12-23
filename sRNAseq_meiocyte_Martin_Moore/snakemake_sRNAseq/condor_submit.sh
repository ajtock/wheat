#!/bin/bash

source activate srna_mapping
snakemake -p --cores 48
source deactivate
