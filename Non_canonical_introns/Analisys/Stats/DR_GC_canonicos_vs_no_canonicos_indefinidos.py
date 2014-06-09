import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC


Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[str(chrfa.id)] = chrfa.seq

	f.close()  

	print >> sys.stderr, "OK"

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = ' ')
	
#	derivados_DR = defaultdict(int)
#	indef_DR = defaultdict(int)
#	canonicos_DR = defaultdict(int)
	
	derivados = []
	indef = []
	canonicos = []			

	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)

#	Canonicos = []
	
#	U2_GTAG = []
#	U12_ATAC = []
#	U12_GTAG = []
#	indef = []

	for row in reader1:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = float(row[8])
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		bodymap_seq = row[16]
		gm12878_seq = row[17]
		DR = int(row[-1])
		
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				derivados.append(row)
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:		
				derivados.append(row)
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				derivados.append(row)		
					
			else:
				indef.append(row)
		
		else:
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				canonicos.append(row)
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				canonicos.append(row)	
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				canonicos.append(row)
				
	print "Indefinidos"

	for row in indef:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])			
		DR = int(row[-1])
		intron_seq = Genome[chr][istart:iend]
		
		print DR, GC(intron_seq)
					
	print "Derivados"
	
	for row in derivados:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])			
		DR = int(row[-1])
		intron_seq = Genome[chr][istart:iend]
		
		print DR, GC(intron_seq)
		
	print "Canonicos"
		
	for row in canonicos:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])			
		DR = int(row[-1])
		intron_seq = Genome[chr][istart:iend]
		
		print DR, GC(intron_seq)

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])	
