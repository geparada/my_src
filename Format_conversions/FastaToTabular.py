import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq



def RandomReads(fasta):

        f = open(fasta)
	SeqTable = []
	table10 = []
	
        for read in SeqIO.parse(f, "fasta"):

		print str(read.id), str(read.seq)


if __name__ == '__main__':
	RandomReads(sys.argv[1])
