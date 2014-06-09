import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

score_U2_GTAG = float(63.2891027914)
score_U12_ATAC = float(60.9280810964)
score_U12_GTAG = float(61.4553595446)

Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def uniform(seq):
	return str(seq).upper()

def main(Final_table):

	for row in csv.reader(open(Final_table), delimiter = '\t'):

		gene = row[0]
		
		intron = row[1]
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		ilength = int(row[6])

		dn = row[7]
		dn_type = row[8]
		dn_type_score = row[9]

		shift = row[26].split(",")
		non_canonical_shift = row[27].split(",")

		L = 15

		exon5 = uniform(Genome[chr][istart-L:istart])
		intron5 = uniform(Genome[chr][istart:istart+L])
		intron3 = uniform(Genome[chr][iend-L:iend])
		exon3 = uniform(Genome[chr][iend:iend+L])		
			
		if strand == "-":
			exon3 = uniform(Genome[chr][istart-L:istart].reverse_complement())
			intron3 = uniform(Genome[chr][istart:istart+L].reverse_complement())
			intron5 = uniform(Genome[chr][iend-L:iend].reverse_complement())
			exon5 = uniform(Genome[chr][iend:iend+L].reverse_complement())


		if shift != ["NO"]: 

			poly5 = 0
			poly3 = 0

			poly_max = 0
			nt_poly_max = ''

			while exon5[-1] == exon5[-1 - poly5]:
				poly5 += 1

			while exon3[0] == exon3[poly3]:
				poly3 += 1

			if exon5[-1] == exon3[0]:
				poly_max = poly5 + poly3
				nt_poly_max = exon5[-1]

			else:
				poly_max, nt_poly_max = max((poly5, exon5[-1]), (poly3, exon3[0]), key=lambda x:x[0])

				
			print gene, intron, exon5, intron5, intron3, exon3, exon5[-1], poly5, exon3[0], poly3, poly_max, nt_poly_max










if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2])	