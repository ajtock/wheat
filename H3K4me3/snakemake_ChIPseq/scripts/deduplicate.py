#!/usr/bin/python2.7
# yupf05@gmail.com
# modified Sebastian Mueller sm934@cam.ac.uk
# call:
# deduplicate.py test_1.fq test_2.fq
# output:
# test_1.RmDup.fq test_2.RmDup.fq

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
import os,sys, argparse


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
	fileName1, fileExtension1 = os.path.splitext(args.input1)
	fileName2, fileExtension2 = os.path.splitext(args.input2)
	outfile1 = open(fileName1+".RmDup"+fileExtension1,"w")
	outfile2 = open(fileName2+".RmDup"+fileExtension2,"w")
	fastq_iter1 = SeqIO.parse(open(args.input1),"fastq")
	fastq_iter2 = SeqIO.parse(open(args.input2),"fastq")
	for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
		if str((rec1+rec2).seq) not in Unique_seqs:
			SeqIO.write(rec1,outfile1,"fastq")
			SeqIO.write(rec2,outfile2,"fastq")
			Unique_seqs.add(str((rec1+rec2).seq))
	outfile1.close()
	outfile2.close()

if __name__ == '__main__':
	Main()
