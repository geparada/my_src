import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def main(bam):

	samfile = pysam.Samfile(bam, "r")
	


	for read in samfile.fetch():
		

		seq = read.seq
		q = read.qual
		if read.flag == 0:			
		
			print "@" + read.qname
			print seq
			print "+"
			print q
		
		else:
			
			print "@" + read.qname
			print Seq(seq).reverse_complement()
			print "+"
			print q[::-1]
			
	samfile.close()



if __name__ == '__main__':
	main(sys.argv[1])

