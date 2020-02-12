#!/bin/bash

/home/ajt200/anaconda3/bin/tabix --preset vcf \
                                 --csi \
                                 -m 14 \
                                 all.GP08_mm75_het3_publication01142019.vcf.gz
