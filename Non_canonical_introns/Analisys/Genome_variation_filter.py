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
	

def DR_generator(intron, ichr, strand, istart, iend, Genome):
	
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
	
	return introns_finded_DR
	

def main(SNP, DVG, Final_table):
	
	reader1 = csv.reader(open(SNP), delimiter = '\t')
	reader2 = csv.reader(open(DVG), delimiter = '\t')
	reader3 = csv.reader(open(Final_table), delimiter = ' ')
	
	variations = defaultdict(set)
	deletion = set([])
	
	for row in reader1:
		chr = row[1]
		start = int(row[2])
		end = int(row[3])
		ID = row[4]
		strand = row[6]
		seq = row[7]
		type = row[11]
		name = chr + ":" + str(start) + strand + str(end) 
		if type=="deletion" or type=="in-del":
			if (end-start)>=6:
				deletion.add((name, chr, strand, start, end, ID))

		
		if start == end:
			variations[chr].add(start)	
		
		else:
			for n in range(start, end):
				variations[chr].add(n)
	
	
	for row in reader2:
		chr = row[1]
		start = row[2]
		end = row[3]
		ID = row[4]
		strand = row[6]
		type = row[10]
		name = chr + ":" + str(start) + strand + str(end) 
		if type=="InDel":
			if (end-start)>=6:
				deletion.add((name, chr, strand, start, end, ID))			
	
	for row in reader3:
	
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = row[8]
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		
		var = variations[chr]
		
		info = " ".join(row)

		introns = set([])

		if bodymap_coverage>=3 and (hg19_cDNA_coverage + hg19_EST_coverage)>=3: 

			confident_intron = True

			for n in range (istart-4, istart+6):
				if (n in var) == True:
					confident_intron = False

			
			for n in range (istart-6, istart+4):
				if (n in var) == True:
					confident_intron = False

			if confident_intron == True:
				introns.add(info)
			
		
		if gm12878_coverage>=3 and (hg19_cDNA_coverage + hg19_EST_coverage)>=3: 
			introns.add(info)
		
		if bodymap_coverage>=3 and gm12878_coverage>=3: 		
			introns.add(info)
		
		if (bodymap_coverage>=3 and (mm9_cDNA_coverage + mm9_EST_coverage)>=3)  or (gm12878_coverage>=3 and (mm9_cDNA_coverage + mm9_EST_coverage)>=3):
			introns.add(info)
		
		deletion_DR = {}
		
		for d in deletion:
			name = row[0]
			chr = row[1]
			strand = row[2]
			start = row[3]
			end = row[4]
			ID = row[5]
			DR= DR_generator(name, chr, strand, int(start), int(end), Genome)
			for i in DR:
				deletion_DR[i] = d
		
		for i in introns:
			#intron = i[0]
			#if deletion_DR.has_key(intron) == False:
			print i

#python ~/my_src/Analisys/Genome_variation_filter.py ~/db/genome/hg19.fa ~/db/Variation/snp135.txt ~/db/Variation/DVG TOTAL_introns.final_table.tags

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2],sys.argv[3],sys.argv[4])
