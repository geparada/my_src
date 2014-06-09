import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def main(fastq):
	for record in SeqIO.parse(fastq, "fastq"):
		print ">" + record.id
		print record.seq
		
if __name__ == '__main__':
	main(sys.argv[1]) 
