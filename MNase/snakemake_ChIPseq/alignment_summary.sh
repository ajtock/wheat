#!/bin/bash

# Usage:
# ./alignment_summary.sh MNase_Rep1 

i=$1

# Specify paths to directories counting read and alignment data
rawDir=./data
dedupDir=./data/dedup
trimmedDir=./data/dedup/trimmed
bamDir=./mapped
bamBothDir=./mapped/both
bamUniqDir=./mapped/unique
statsDir=./logs/alignment_stats

[ -d ./logs/alignment_stats ] || mkdir ./logs/alignment_stats

zcat $rawDir/${i}_R1.fastq.gz | echo $((`wc -l`/4)) > $statsDir/${i}_R1.fastq.gz.stats
zcat $rawDir/${i}_R2.fastq.gz | echo $((`wc -l`/4)) > $statsDir/${i}_R2.fastq.gz.stats
zcat $dedupDir/${i}_R1_dedup_repair.fastq.gz | echo $((`wc -l`/4)) > $statsDir/${i}_R1_dedup_repair.fastq.gz.stats
zcat $dedupDir/${i}_R2_dedup_repair.fastq.gz | echo $((`wc -l`/4)) > $statsDir/${i}_R2_dedup_repair.fastq.gz.stats
zcat $trimmedDir/${i}_R1_dedup_repair_trimmed.fastq.gz | echo $((`wc -l`/4)) > $statsDir/${i}_R1_dedup_repair_trimmed.fastq.gz.stats
zcat $trimmedDir/${i}_R2_dedup_repair_trimmed.fastq.gz | echo $((`wc -l`/4)) > $statsDir/${i}_R2_dedup_repair_trimmed.fastq.gz.stats

samtools view $bamDir/${i}_MappedOn_wheat_v1.0.bam | cut -f1 | sort | uniq | wc -l > $statsDir/${i}_MappedOn_wheat_v1.0.bam.stats
samtools view $bamBothDir/${i}_MappedOn_wheat_v1.0_lowXM_both_sort.bam | cut -f1 | sort | uniq | wc -l > $statsDir/${i}_MappedOn_wheat_v1.0_lowXM_both_sort.bam.stats
samtools view $bamUniqDir/${i}_MappedOn_wheat_v1.0_lowXM_unique_sort.bam | cut -f1 | sort | uniq | wc -l > $statsDir/${i}_MappedOn_wheat_v1.0_lowXM_unique_sort.bam.stats

paste -d "\t" $statsDir/${i}_R1.fastq.gz.stats $statsDir/${i}_R1_dedup_repair.fastq.gz.stats $statsDir/${i}_R1_dedup_repair_trimmed.fastq.gz.stats $statsDir/${i}_MappedOn_wheat_v1.0.bam.stats $statsDir/${i}_MappedOn_wheat_v1.0_lowXM_both_sort.bam.stats $statsDir/${i}_MappedOn_wheat_v1.0_lowXM_unique_sort.bam.stats > $statsDir/${i}_alignment_summary.txt
sed -i '1i Total_sequenced_read_pairs\tDeduplicated\tTrimmed\tPrimary_alignments\tMAPQ>=2_mismatches<=6\tMAPQ>=23_mismatches<=6' $statsDir/${i}_alignment_summary.txt
