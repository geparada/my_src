import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def FastaToTabular(fasta):
        f = open(fasta)
	SeqTable = []

        for chrfa in SeqIO.parse(f, "fasta"):
                print str(chrfa.id) + "\t" +  str(chrfa.seq)
		

	f.close()


if __name__ == '__main__':
        FastaToTabular(sys.argv[1])

