import sys
import csv
import re
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

def main(TOTAL_Tissues_introns):

	for intron in csv.reader(open(TOTAL_Tissues_introns)):
		chr, istart, iend = re.findall(r"[\w']+", intron[0])

		istart = int(istart)
		iend = int(iend)

		dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]

		if '-' in intron[0]:
			dn = dn.reverse_complement()

		dn = str(dn).upper()

		if dn!='GTAG' and dn!='GCAG' and dn!='ATAC':

			print dn

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2])	