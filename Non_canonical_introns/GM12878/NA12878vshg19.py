import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable=[]
csv.field_size_limit(1000000000)

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

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
	
	return introns_finded_DR
	
	
def main(NA12878, hg19):
	reader1 = csv.reader(open(NA12878), delimiter = ' ')
	reader2 = csv.reader(open(hg19), delimiter = ' ')
	
	intron_list = set([])
	Genome = dict(SeqTable)	
	
	NA12878_list = []
	hg19_list = []
	
	
	for row in reader1:
		intron = row[0]
		coverage = row[1]
		chr = row[2]
		strand = row[3]
		start = row[4]
		end = row[5]
		dn = row[7]
		DR_introns = DR_counter(intron, chr, strand, int(start), int(end), Genome)
		DR_intron_ID = ','.join(sorted(DR_introns))
		#if dn != 'GTAG' and dn != 'GCAG' and dn != 'ATAC':
		NA12878_list.append((DR_intron_ID, [intron, chr, strand, start, end, dn, coverage]))
		intron_list.add((DR_intron_ID))
		
	NA12878_dict = dict(NA12878_list)
		
	for row in reader2:
		hg19_intron = row[0]
		hg19_coverage = row[1]
		hg19_chr = row[2]
		hg19_strand = row[3]
		hg19_start = row[4]
		hg19_end = row[5]
		hg19_dn = row[7]
		DR_introns = DR_counter(hg19_intron, hg19_chr, hg19_strand, int(hg19_start), int(hg19_end), Genome)
		DR_intron_ID = ','.join(sorted(DR_introns))
		if hg19_dn != 'GTAG' and hg19_dn != 'GCAG' and hg19_dn != 'ATAC':
		#	hg19_list.append((DR_intron_ID, [intron, chr, strand, start, end, dn, coverage]))
		#	intron_list.add((DR_intron_ID))

			try:
				info = NA12878_dict[DR_intron_ID]

				#if bodymap_dict.has_key(DR_intron_ID)== True:
				NA12878_intron = info[0]
				NA12878_chr = info[1]
				NA12878_strand = info [2]
				NA12878_start = info[3]
				NA12878_end = info[4]
				NA12878_dn = info[5]		
				NA12878_coverage = info[6]
				if NA12878_dn == 'GTAG' or NA12878_dn == 'GCAG' or NA12878_dn == 'ATAC':
					print hg19_intron, "Canonical_in_hg19", hg19_dn, NA12878_dn
				
					
			except KeyError:
				print hg19_intron, "NA12878_only", hg19_dn, "NO"
				#pass	





		
	#NA12878_dict = dict(NA12878_list)
	#hg19_dict = dict(hg19_list)
	
	#intron_list_sort = sorted(list(intron_list)) 
	
	
	
	

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2],sys.argv[3])
