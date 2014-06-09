from random import random, choice
import sys
from itertools import groupby
import csv

def fasta_iter(fasta_name):
	"""
	given a fasta file. yield tuples of header, sequence
	"""
	fh = open(fasta_name)
	# ditch the boolean (x[0]) and just keep the header or sequence since
	# we know they alternate.
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
	for header in faiter:
		# drop the ">"
		header = header.next()[1:].strip()
		# join all sequence lines to one.
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq

def main(fasta_name, mutation_freq_index):

	reader = csv.reader(open(mutation_freq_index), delimiter = ' ')
	nt_P = []	
	for row in reader:
		nt = row[0]      #Nucleotide position number          
		P = row[2]       # Error probability
		nt_P.append((nt, P))

	index = dict(nt_P)
#	for i in range(10):
#		print index[str(i)] 

	for header, seq in fasta_iter(fasta_name):
		seq = list(seq)
		for i, s in enumerate(seq):
			val = random()
			mutation_freq = float(index[str(i)])
			if val < mutation_freq:
				#choose a random nucleotide that's different.
				seq[i] = choice([x for x in "ACTG" if  x != s.upper()])
		
		print ">%s\n%s" % (header, "".join(seq))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
