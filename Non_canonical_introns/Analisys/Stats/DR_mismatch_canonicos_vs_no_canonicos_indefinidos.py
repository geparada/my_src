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

def DR_mismatch_counter(intron, ichr, strand, istart, iend, Genome):
	
	
	max_mismatchs = 1
			
	L = 100         #Solo permite que se corra L pares de bases para buscar DR 		

	#Extrayendo regiones exonicas colindantes
			
	SJ5U = Genome[ichr][istart-L : istart].lower()
	SJ5D = Genome[ichr][istart : istart+L].lower()
	SJ3U = Genome[ichr][iend-L : iend].lower()
	SJ3D = Genome[ichr][iend : iend+L].lower()
	
	if strand == "-":
		SJ5U = Genome[ichr][iend : iend+L].lower().reverse_complement() 
		SJ5D = Genome[ichr][iend-L : iend].lower().reverse_complement()
		SJ3U = Genome[ichr][istart : istart+L].lower().reverse_complement()
		SJ3D = Genome[ichr][istart-L : istart].lower().reverse_complement()
	
	DRU = 0
	DRD = 0
								
	#Contando directos repetidos y generando intrones no consenso alternativos
	
	mismatchs = 0
	i = 0
		
	try:
		while mismatchs <= max_mismatchs:
			i += 1
			
			if  SJ5U[L-1-i]==SJ3U[L-1-i]: 
				DRU += 1
			else:
				mismatchs += 1
				
			if mismatchs > max_mismatchs:	
				break
				
	except IndexError:
		pass 
	
	mismatchs = 0	
	i = 0
		
	try:
		while mismatchs <= max_mismatchs:
			i += 1
			
			if SJ5D[i]==SJ3D[i]:
				DRD += 1
			else:
				mismatchs += 1
				
			if mismatchs > max_mismatchs:	
				break
	except IndexError:
		pass
	
	return DRU + DRD

def main(Final_table):
	
	derivados_DR = defaultdict(int)
	indef_DR = defaultdict(int)
	canonicos_DR = defaultdict(int)
	
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
		DR = DR_mismatch_counter(intron, chr, strand, int(istart), int(iend), Genome)
				
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				derivados_DR[DR] += 1
				derivados += 1
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				derivados_DR[DR] += 1		
				derivados += 1
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				derivados_DR[DR] += 1
				derivados += 1				
					
			else:
				indef_DR[DR] += 1
				indef += 1
		
		else:
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				canonicos_DR[DR] += 1
				canonicos += 1
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				canonicos_DR[DR] += 1
				canonicos += 1		
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				canonicos_DR[DR] += 1				
				canonicos += 1		
		
	for i in range(201):
		print i, percent(indef_DR[i], indef) , percent(derivados_DR[i], derivados), percent(canonicos_DR[i], canonicos) 


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])
