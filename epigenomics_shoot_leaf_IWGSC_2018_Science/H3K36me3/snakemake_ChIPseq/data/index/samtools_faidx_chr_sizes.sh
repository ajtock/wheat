#!/bin/bash

i=$1

samtools faidx ${i}.fa
cut -f1,2 ${i}.fa.fai > ${i}.fa.sizes
