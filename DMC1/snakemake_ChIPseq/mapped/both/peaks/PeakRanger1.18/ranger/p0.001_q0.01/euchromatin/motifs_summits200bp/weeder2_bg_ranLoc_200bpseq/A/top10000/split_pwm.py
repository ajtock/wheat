#!/home/ajt200/anaconda3/bin/python3

# Parse weeder2-generated sequence of position weight matrices (PWMs)
# and split into individual PWM files

# Usage:
# python3 ./split_pwm.py ./DMC1_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp.fa.matrix.w2

import os
from os import path
import sys

# Set/get current working directory
#os.chdir("/home/ajt200/analysis/wheat/DMC1/snakemake_ChIPseq/mapped/both/peaks/PeakRanger1.18/ranger/p0.001_q0.01/euchromatin/motifs_summits200bp/weeder2_bg_ranLoc_200bpseq/A/top10000")
print(os.getcwd())
#infile = open("./DMC1_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp.fa.matrix.w2")
#matrix_file = "./DMC1_Rep1_ChIP_peaks_minuslog10Qsorted_in_Agenome_euchromatin_summits200bp.fa.matrix.w2"
infile = open(sys.argv[1])
matrix_file = sys.argv[1]

#path = matrix_file.split('/')[0]
path = sys.argv[1].split('/')[0]

print('\n' + "matrix file: " + matrix_file + '\n')
num_pwms = 'grep "^>" %s | wc -l' % (matrix_file)
print("Number of pwms in matrix file = ")
os.system(num_pwms)

# Assume outfile is not open
opened = False

for line_matrix in infile:
    # If line begins with ">"
    if line_matrix[0] == ">":
        # Close the preceding outfile if it is open (see below and follow loop)
        if(opened):
            outfile.close()
        # Set opened to True to represent an opened outfile
        opened = True
        # Extract pwm name: remove ">", extract pwm string (replacing tab with underscore),
        # and remove any spaces or new lines
        # (e.g., ">MAT1   TACATTTTTT" will become "MAT1_TACATTTTTT")
        pwm_name = line_matrix[1:].rstrip().replace("\t", "_")
        print("pwm: " + pwm_name)
        print("Name of outfile: " + path + "/" + str(pwm_name) + ".pwm")
        outfile = open(path + "/" + str(pwm_name) + ".pwm", 'w')
    # Write the line to the file.
    # If the line begins with ">" (at position 0; i.e., is a PWM name),
    # a new outfile is created and opened, ready to be written to.
    # Subsequent lines (i.e., the rows of the PWMs containing base relative frequencies)
    # are each successively written to the outfile with each loop.
    # Only when the for-loop reaches a new line beginning with ">" is the outfile closed,
    # and a new outfile with a new name is created and opened.
    outfile.write(line_matrix)
outfile.close()
print("Finished")
