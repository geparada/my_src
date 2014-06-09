import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

csv.field_size_limit(1000000000)

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
		Genome[str(chrfa.id)] = chrfa.seq

	f.close()  

	print >> sys.stderr, "OK"

def SJ5_3_identity(intron, ichr, strand, istart, iend, Genome):
	
			
	L = 5		

	I = 0
	
	SJ5 = Genome[ichr][istart - L : istart + L].lower()
	SJ3 = Genome[ichr][iend - L : iend + L].lower()
	
	if strand == "-":
		SJ5 = Genome[ichr][iend - L : iend + L].lower().reverse_complement() 
		SJ3 = Genome[ichr][istart - L : istart + L].lower().reverse_complement()
			
	for i in range(len(SJ5)):
		if SJ5[i] == SJ3[i]:
			I += 1
		
	return percent(I, 2*L)

def main(Final_table):
	
	derivados_PI = defaultdict(int)
	indef_PI = defaultdict(int)
	canonicos_PI = defaultdict(int)
	
	derivados = 0
	indef = 0
	canonicos = 0			

	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)
	
	for row in csv.reader(open(Final_table), delimiter = '\t'):
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
		PI = int(str(SJ5_3_identity(intron, chr, strand, int(istart), int(iend), Genome)/100)[2]) + 10*int(str(SJ5_3_identity(intron, chr, strand, int(istart), int(iend), Genome)/100)[0])
 
				
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				derivados_PI[PI] += 1
				derivados += 1
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				derivados_PI[PI] += 1		
				derivados += 1
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				derivados_PI[PI] += 1
				derivados += 1				
					
			else:
				indef_PI[PI] += 1
				indef += 1
		
		else:
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				canonicos_PI[PI] += 1
				canonicos += 1
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				canonicos_PI[PI] += 1
				canonicos += 1		
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				canonicos_PI[PI] += 1				
				canonicos += 1		
		
	for i in range(11):
		print i*10, percent(indef_PI[i], indef) , percent(derivados_PI[i], derivados), percent(canonicos_PI[i], canonicos)
	
	 


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])
