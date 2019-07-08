#!/bin/bash

# Usage:
# ./pdfTotiffTopng.sh ASY1_CS_peaks_in_Agenome_euchromatin ASY1_CS_Rep1_ChIP bodies

features=$1
factor=$2
region=$3

#[ -d tiff ] || mkdir tiff
[ -d png ] || mkdir png

for i in *"_around_"${features}"_heatmap_ordered_by_log2_"${factor}"_control_in_"${region}".pdf"
do
#( gs -q -dNOPAUSE -r300x300 -sDEVICE=tiff24nc -sOutputFile=${i}.tiff ${i} -c quit
(  convert -density 300 ${i} -quality 90 ${i}.png ) &
done
wait

#rename -v 's/.pdf.tiff/.tiff/' *.tiff
rename -v 's/.pdf.png/.png/' *.png

#mv *.tiff tiff
mv *.png png
