import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict


score_U2_GTAG = float(63.2891027914)
score_U12_ATAC = float(60.9280810964)
score_U12_GTAG = float(61.4553595446)

Genome = {}

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

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

	indef_shift = defaultdict(int)
	indef_non_shift = defaultdict(int)
	canonicos = defaultdict(int)

	total_indef_shift = 0
	total_indef_non_shift = 0
	total_canonicos = 0


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
		dn_type_score = float(row[9])

		shift = row[26].split(",")
		non_canonical_shift = row[27].split(",")

		L = 50

		exon5 = uniform(Genome[chr][istart-L:istart])
		intron5 = uniform(Genome[chr][istart:istart+L])
		intron3 = uniform(Genome[chr][iend-L:iend])
		exon3 = uniform(Genome[chr][iend:iend+L])		
			
		if strand == "-":
			exon3 = uniform(Genome[chr][istart-L:istart].reverse_complement())
			intron3 = uniform(Genome[chr][istart:istart+L].reverse_complement())
			intron5 = uniform(Genome[chr][iend-L:iend].reverse_complement())
			exon5 = uniform(Genome[chr][iend:iend+L].reverse_complement())

		poly5 = 0
		poly3 = 0
		poly_max = 0
		nt_poly_max = ''

		while exon5[-1] == exon5[-1 - poly5]:
			poly5 += 1
			if poly5 == L:
				break
		while exon3[0] == exon3[poly3]:
			poly3 += 1
			if poly3 == L:
				break

		if exon5[-1] == exon3[0]:
			poly_max = poly5 + poly3
			nt_poly_max = exon5[-1]
		else:
			poly_max, nt_poly_max = max((poly5, exon5[-1]), (poly3, exon3[0]), key=lambda x:x[0])


		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				pass
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:		
				pass
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				pass	
					
			else:

				if shift != ["NO"]:
					#print "\t".join(row)
					if poly_max >= 5:
						print row

					indef_shift[poly_max] += 1
					total_indef_shift += 1

				else:
					#print "\t".join(row)

					indef_non_shift[poly_max] += 1
					total_indef_non_shift += 1

		else:
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				canonicos[poly_max] += 1
				total_canonicos += 1
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				canonicos[poly_max] += 1
				total_canonicos += 1
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				canonicos[poly_max] += 1

#	print "indef_shift", "indef_non_shift", "canonicos"

#	for i in range(1,L+1):
#		print i, percent(indef_shift[i], total_indef_shift) , percent(indef_non_shift[i], total_indef_non_shift), percent(canonicos[i], total_canonicos) 

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2])	
