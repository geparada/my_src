import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

hexamers = defaultdict(int)


def main (fasta):
	""" Realiza una tabla de frecuencias de los k-meks precentes """

	for record in SeqIO.parse(fasta, "fasta"):
		kmer_len = 6
		nhexamers = 0
		
		for i in range(len(record.seq)- kmer_len):
			hexamers[str(record.seq[i:i+kmer_len])] += 1
			nhexamers += 1
			
	hexamers_list = hexamers.items()

	for k in sorted(hexamers_list, key=lambda a: a[1], reverse=True):
		print k[0], k[1], float(k[1])/float(nhexamers) * 100



if __name__ == '__main__':
	main(sys.argv[1])	



