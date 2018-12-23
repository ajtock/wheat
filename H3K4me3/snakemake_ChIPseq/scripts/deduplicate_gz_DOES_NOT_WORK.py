#!/usr/bin/python2.7
# yupf05@gmail.com
# modified by Sebastian Mueller (sm934@cam.ac.uk)
# modified by Andy Tock (ajt200@cam.ac.uk) for use with fastq.gz files
### ajt changes based on https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz
### and https://bioinformatics.stackexchange.com/questions/892/how-do-you-write-a-gz-fastq-file-with-biopython/894
# call:
# deduplicate_gz.py test_1.fastq.gz test_2.fastq.gz
# output:
# test_1.RmDup.fastq.gz test_2.RmDup.fastq.gz

from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
import itertools
import os,sys, argparse
import gzip

def ParseArg():
	p=argparse.ArgumentParser( description = 'Remove duplicated reads which have same sequences for both forward and reverse reads. Choose the one appears first.', epilog = 'Library dependency: Bio, itertools')
	p.add_argument('input1',type=str,metavar='reads1',help='forward input fastq/fasta file')
	p.add_argument('input2',type=str,metavar='reads2',help='reverse input fastq/fasta file')
	if len(sys.argv)==1:
		print >>sys.stderr,p.print_help()
		exit(0)
	return p.parse_args()


def Main():

	Unique_seqs=set()
	args=ParseArg()
        fileName1 = os.path.basename(args.input1)
        fileName2 = os.path.basename(args.input2)
        indexOfDot1 = fileName1.index(".")
        indexOfDot2 = fileName2.index(".")
        filePrefix1 = fileName1[:indexOfDot1]
        filePrefix2 = fileName2[:indexOfDot2]
	outfile1 = gzip.open(filePrefix1+".RmDup.fastq.gz","wb")
	outfile2 = gzip.open(filePrefix2+".RmDup.fastq.gz","wb")
	fastq_iter1 = SeqIO.parse(gzip.open(args.input1, "rb"),"fastq")
	fastq_iter2 = SeqIO.parse(gzip.open(args.input2, "rb"),"fastq")
	for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
		if str((rec1+rec2).seq) not in Unique_seqs:
			SeqIO.write(rec1,outfile1,"fastq")
			SeqIO.write(rec2,outfile2,"fastq")
			Unique_seqs.add(str((rec1+rec2).seq))
	outfile1.close()
	outfile2.close()

if __name__ == '__main__':
	Main()
