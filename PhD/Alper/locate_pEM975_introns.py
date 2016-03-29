import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re


def main(pEM975, pEM975_introns):

	f1 = open(pEM975)
	f2 = open(pEM975_introns)

	for pEM975_fa in SeqIO.parse(f1, "fasta"):


		c = 0
		estart = 0

		for introns_fa in SeqIO.parse(f2, "fasta"):


			istart = re.search(str(introns_fa.seq), str(pEM975_fa.seq)).start()
			iend = istart + len(introns_fa.seq)


			eend = istart

			#print "exon", estart, eend, pEM975_fa.seq[estart:eend]
			#print "intron", istart, iend, pEM975_fa.seq[istart:iend]

			#out = ["pEM975", "CE10_exonic", "exon", estart+1, eend+1, ".", "+", ".", "gene_id pEM975_gene;"]
			out = ["pEM975", "CE10_exonic", "exon", istart+1, iend+1, ".", "+", ".", "gene_id pEM975_gene;"]

			estart = iend

			

			print "\t".join(map(str, out))






if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])