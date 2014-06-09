import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()
	
def main (chr, start, end, strand):
	""" Extrae genoma de coordenadas """
	
	ID = chr + ":" + str(start) + strand + str(end)
	seq = Genome[chr][start:end]
	if strand == "-":
		seq = seq.reverse_complement()

	print ">" + ID
	print str(seq)


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])	
