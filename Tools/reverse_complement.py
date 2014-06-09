import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def main (seq):
	print seq.reverse_complement()
			


if __name__ == '__main__':
	main(Seq(sys.argv[1]))
