#!/bin/bash

# Usage:
# csmit -m 5G -c 1 "bash ./plotFingerprint.sh H2AZ_Rep1_ChIP H3K4me3_Rep1_ChIP H3K4me3_ChIP_SRR6350668 H3K9me2_Rep1_ChIP H3K9ac_ChIP_SRR6350667 H3K27me1_Rep1_ChIP H3K27me3_ChIP_SRR6350666 H3K36me3_ChIP_SRR6350670 CENH3_ChIP_SRR1686799 MNase_Rep1 H3_input_SRR6350669 max"

dat1=${1}
dat2=${2}
dat3=${3}
dat4=${4}
dat5=${5}
dat6=${6}
dat7=${7}
dat8=${8}
dat9=${9}
dat10=${10}
dat11=${11}
threads=${12}

plotFingerprint -b bam/${dat1}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat2}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat3}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat4}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat5}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat6}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat7}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat8}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat9}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat10}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                   bam/${dat11}_MappedOn_wheat_v1.0_lowXM_both_sort.bam \
                -o wheat_ChIP_input_fingerprints.pdf \
                --extendReads \
                --outRawCounts wheat_ChIP_input_fingerprints.tab \
                -p ${threads}
