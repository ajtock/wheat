#!/bin/bash

# Usage
# ./awk_download_list.sh PRJNA420988

accession=$1

awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' $accession".txt" \
| cut -f 3 \
| awk -F ";" 'OFS="\n" {print $1, $2}' \
| awk NF \
| awk 'NR > 1, OFS="\n" { print "ascp -QT -l 300m -P33001 -i /home/ajt200/.aspera/connect/etc/asperaweb_id_dsa.openssh" " " $1 " ."}' \
> $accession"_download.txt"
