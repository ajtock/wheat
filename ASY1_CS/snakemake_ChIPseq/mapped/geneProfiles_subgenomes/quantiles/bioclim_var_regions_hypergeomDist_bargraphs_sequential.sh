#!/bin/bash

for i in {01..19}
do
( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
    'HudsonRM_all' 'bodies' 1 2 'genomewide' \
     BIO${i}_LAR_overlapping BIO${i}' LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
done
./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
  'HudsonRM_all' 'bodies' 1 2 'genomewide' \
  Altitude_LAR_overlapping 'Altitude LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen'

#for i in {01..19}
#do
#( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#    'HudsonRM_all' 'bodies' 1 2 'genomewide' \
#     BIO${i}_LAR2_overlapping BIO${i}' LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
#done
#./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#  'HudsonRM_all' 'bodies' 1 2 'genomewide' \
#  Altitude_LAR2_overlapping 'Altitude LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
#
#for i in {01..19}
#do
#( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#    'ASY1_CS_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#     BIO${i}_LAR_overlapping BIO${i}' LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
#done
#./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#  'ASY1_CS_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#  Altitude_LAR_overlapping 'Altitude LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
#
#for i in {01..19}
#do
#( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#    'ASY1_CS_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#     BIO${i}_LAR2_overlapping BIO${i}' LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
#done
#./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#  'ASY1_CS_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#  Altitude_LAR2_overlapping 'Altitude LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
#
#for i in {01..19}
#do
#( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#    'DMC1_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#     BIO${i}_LAR_overlapping BIO${i}' LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
#done
#./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#  'DMC1_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#  Altitude_LAR_overlapping 'Altitude LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
#
#for i in {01..19}
#do
#( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#    'DMC1_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#     BIO${i}_LAR2_overlapping BIO${i}' LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
#done
#./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#  'DMC1_Rep1_ChIP' 'bodies' 1 4 'genomewide' \
#  Altitude_LAR2_overlapping 'Altitude LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
#
#for i in {01..19}
#do
#( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#    'cMMb' 'bodies' 1 4 'genomewide' \
#     BIO${i}_LAR_overlapping BIO${i}' LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
#done
#./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#  'cMMb' 'bodies' 1 4 'genomewide' \
#  Altitude_LAR_overlapping 'Altitude LAR-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
#
#for i in {01..19}
#do
#( ./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#    'cMMb' 'bodies' 1 4 'genomewide' \
#     BIO${i}_LAR2_overlapping BIO${i}' LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen' ) &
#done
#./proportion_query_genes_in_gene_quantiles_hypergeometricTest_bargraph_only.R \
#  'cMMb' 'bodies' 1 4 'genomewide' \
#  Altitude_LAR2_overlapping 'Altitude LAR2-overlapping' 100000 'black,darkgreen,seagreen,springgreen'
