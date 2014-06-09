from random import random, choice
import sys
from itertools import groupby
import csv
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


def main(fastq_name):

	MAX_Q = 100
	logQs = [10 ** (float(-Q) / 10) for Q in range(MAX_Q + 1)]
	for record in SeqIO.parse(fastq_name, "fastq"):
	
		#print _get_phred_quality(record)
			
	#for header, seq in fasta_iter(fasta_name):

		seq = list(record.seq)
		Q = record.letter_annotations["phred_quality"]
		P = [10 ** (float(-x) / 10) for x in Q]
		
		
	

		for i, s in enumerate(seq):
			val = random()
			mutation_freq = P[i]
			if val < mutation_freq:
				#choose a random nucleotide that's different.
				seq[i] = choice([x for x in "ACTG" if  x != s.upper()])

		#record.seq = "".join(seq)
		#print record.format("fastq")
		
		#print "@" + str(record.ide) + "\n" + "".join(seq) + "\n" + "+" + "\n" + Q
		mutated = SeqRecord( Seq("".join(seq), generic_dna), id = record.id, description = "" )
		mutated.letter_annotations["phred_quality"] = Q
		print mutated.format("fastq"),
		 

if __name__ == "__main__":
    main(sys.argv[1])
