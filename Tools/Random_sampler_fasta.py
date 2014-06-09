import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import random


def main (fasta, N):
	""" Seleciona un numero random de secuencias"""
	
	FASTA = []
	
	for record in SeqIO.parse(fasta, "fasta"):
		
		FASTA.append((record.id, str(record.seq)))
	
	for record in random.sample(FASTA, N):
		id = record[0]
		seq = record[1]
		
		print ">" + id
		print seq



if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2])) 
