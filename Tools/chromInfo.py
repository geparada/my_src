import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from random import randint, sample
from operator import itemgetter
from collections import defaultdict
from operator import itemgetter

chrom_sizes = defaultdict(int)

def main(genome_fasta):

	for record in SeqIO.parse(genome_fasta, "fasta"):

		chrom_sizes[record.id] += len(record.seq)


	for i in chrom_sizes.items():

		print "\t".join(map(str, i))


if __name__ == '__main__':
	main(sys.argv[1])