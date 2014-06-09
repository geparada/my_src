import sys
import csv
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


def DR_counter(intron, ichr, strand, istart, iend, Genome):
	
	introns_finded_DR = [intron]
			
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
		
	try:
		while SJ5U[L-1-DRU]==SJ3U[L-1-DRU]:
			DRU += 1
			if strand == "+":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart-DRU) + strand + str(iend-DRU)]
			elif strand == "-":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart+DRU) + strand + str(iend+DRU)]
			if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
				break
	except IndexError:
		pass 
	try:
		while SJ5D[DRD]==SJ3D[DRD]:
			DRD += 1
			if strand == "+":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart+DRD) + strand + str(iend+DRD)]
			elif strand == "-":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart-DRD) + strand + str(iend-DRD)]
			if SJ5D[DRD]!=SJ3D[DRD]:
				break
	except IndexError:
		pass
	
	return (DRU, DRD)

def main(SNP, Final_table):
	
	csv.field_size_limit(sys.maxsize)
	
	reader1 = csv.reader(open(SNP), delimiter = '\t')
	reader2 = csv.reader(open(Final_table), delimiter = ' ')
	
	SNPs = defaultdict(set)
	
	for row in reader1:
		chr = row[1]
		start = int(row[2])
		end = int(row[3])
		ID = row[4]
		strand = row[6]
		seq = row[7]
		type = row[11]
		name = chr + ":" + str(start) + strand + str(end) 
		
		if start == end:
			SNPs[chr].add(start)	
		
		else:
			for n in range(start, end):
				SNPs[chr].add(n)	
	
	for row in reader2:
	
		intron = row[0]
		coverage = int(row[1])
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		ilength = row[6]
		dn = row[7]
		dn_type = row[8]
		dn_type_score = row[9]

		
		chr_SNPs = SNPs[chr]

		introns = set([])

		confident_intron = True
		DR_intron = DR_counter(intron, chr, strand, int(istart), int(iend), Genome)
		DRU = DR_intron[0]
		DRD = DR_intron[1]
		
		#if bodymap_coverage>=3 and (hg19_cDNA_coverage + hg19_EST_coverage)>=3 and gm12878_coverage < 3 and (mm9_cDNA_coverage + mm9_EST_coverage)<3: 

		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":

			if strand == "+":

				for n in range (istart-1-DRU, istart+3+DRD):
					if (n in chr_SNPs) == True:
						confident_intron = False
						
				for n in range (iend-3-DRU, iend+1+DRD):
					if (n in chr_SNPs) == True:
						confident_intron = False
							
			if strand == "-":
					
				for n in range (istart-1-DRD, istart+3+DRU):
					if (n in chr_SNPs) == True:
						confident_intron = False

				for n in range (iend-3-DRD, iend+1+DRU):
					if (n in chr_SNPs) == True:
						confident_intron = False
		
		if confident_intron:
			print " ".join(row)
		

		

#python ~/my_src/cDNA_EST/SNPs_filter.py ~/db/genome/hg19.fa ~/db/Variation/snp135.txt non_canonical.final_table.tags.seq.DR

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2],sys.argv[3])
