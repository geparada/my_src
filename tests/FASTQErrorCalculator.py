import sys
from Bio import SeqIO
from Bio.Seq import Seq




def ErrorCal(fastqfile):
	
	f = open(fastqfile)
	Psum = 0
	Nsum = 0

	
	for record in SeqIO.parse(f, "fastq"):
		Nsum = Nsum + len(record.seq)		

		for Q in record.letter_annotations["phred_quality"]:
			P = 10**(float(-Q)/10)
			Psum = Psum + P 

	print Psum/Nsum
		



if __name__ == '__main__':
	ErrorCal(sys.argv[1])   
