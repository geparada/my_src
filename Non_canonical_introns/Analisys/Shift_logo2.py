import sys
import csv
from collections import defaultdict
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[str(chrfa.id)] = chrfa.seq

	f.close()  

	print >> sys.stderr, "OK"


def main (FINAL_TABLE):
	""" Extrae secuencias par logo de los shift """
	
	reader1 = csv.reader(open(FINAL_TABLE), delimiter="\t")
	
	for row in reader1:
		
		intron = row[1]
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		shifts = row[26].split(",")
		
		shift_1_2 = "NO"
		shift_2_1 = "NO"
		
		for i in shifts:
				alt_intron, dn, fold, donor, aceptor = i.split("|")
				donor = int(donor)
				aceptor = int(aceptor)
				
				if donor == 1 and aceptor == 2:
					shift_1_2 = "YES"
				
				if donor == -2 and aceptor == -1:
					shift_2_1 = "YES"
		
		if shift_2_1 == "YES":
			print Genome[chr][istart-15:istart+15], Genome[chr][iend-15:iend+15]
#			  Genome[chr][(iend-2):iend]		

					
	
	



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2])
