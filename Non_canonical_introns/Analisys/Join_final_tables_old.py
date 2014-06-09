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
	

def main(bodymap, GM17828, genecode, cDNA_EST):
	reader1 = csv.reader(open(bodymap), delimiter = ' ')
	reader2 = csv.reader(open(GM17828), delimiter = ' ')
	reader3 = csv.reader(open(genecode), delimiter = ' ')
	reader4 = csv.reader(open(cDNA_EST), delimiter = '\t')
	
	intron_list = set([])
	Genome = dict(SeqTable)	
	
	bodymap_list = []
	GM17828_list = []
	genecode_list = []
	cDNA_EST_list = []
	
	
	
	for row in reader1:
		intron = row[0]
		coverage = row[1]
		chr = row[2]
		strand = row[3]
		start = row[4]
		end = row[5]
		dn = row[7]
		DR_introns = DR_counter(intron, chr, strand, int(start), int(end), Genome)
		DR_intron_ID = ','.join(sorted(DR_introns))               #Esta es la manera de que dos intrones no canonicos iguales siempre tengan el mismo ID
		if dn != 'GTAG' and dn != 'GCAG' and dn != 'ATAC':
			bodymap_list.append((DR_intron_ID, [intron, chr, strand, start, end, dn, coverage]))
			intron_list.add((DR_intron_ID))
		
	for row in reader2:
		intron = row[0]
		coverage = row[1]
		chr = row[2]
		strand = row[3]
		start = row[4]
		end = row[5]
		dn = row[7]
		DR_introns = DR_counter(intron, chr, strand, int(start), int(end), Genome)
		DR_intron_ID = ','.join(sorted(DR_introns))
		if dn != 'GTAG' and dn != 'GCAG' and dn != 'ATAC':
			GM17828_list.append((DR_intron_ID, [intron, chr, strand, start, end, dn, coverage]))
			intron_list.add((DR_intron_ID))
		
	for row in reader3:
		transcript = row[0]
		chr = row[1]
		start = row[2]
		end = row[3]
		strand = row[4]
		intron = row[6]
		dn = row[7]
		DR_introns = DR_counter(intron, chr, strand, int(start), int(end), Genome)
		DR_intron_ID = ','.join(sorted(DR_introns))
		if dn != 'GTAG' and dn != 'GCAG' and dn != 'ATAC':
			genecode_list.append((DR_intron_ID, [intron,chr, strand, start, end, dn,transcript]))
			intron_list.add((DR_intron_ID))
			
	for row in reader4:
		dn = row[1]
		chr = row[2]
		start = row[3]
		end = row[4]
		strand = row[5]
		intron = chr + ":" + start + strand + end
		coverage_EST = row[7]
		coverage_cDNA = row[8]
		EST = row[9]
		cDNA = row[10]
		DR_introns = DR_counter(intron, chr, strand, int(start), int(end), Genome)
		DR_intron_ID = ','.join(sorted(DR_introns))
		if dn != 'GTAG' and dn != 'GCAG' and dn != 'ATAC':
			cDNA_EST_list.append((DR_intron_ID, [intron,chr, strand, start, end, dn, coverage_EST,coverage_cDNA,EST,cDNA]))
			intron_list.add((DR_intron_ID))
	
	
	
	bodymap_dict = dict(bodymap_list)
	GM17828_dict = dict(GM17828_list)
	genecode_dict = dict(genecode_list)
	cDNA_EST_dict = dict(cDNA_EST_list)		
			
	#keys = 	set(intron_list[0])

#	print intron_list

	intron_list_sort = sorted(list(intron_list))              #Este paso solo sirve para que queden más ordenados.
			
	for row in intron_list_sort:
		DR_intron_ID = row
		
		chr = ''
		strand = ''
		start = ''
		end = ''
		dn = ''
	
		try:
			info = cDNA_EST_dict[DR_intron_ID]
			
			if cDNA_EST_dict.has_key(DR_intron_ID)== True:
				intron = info[0]
				chr = info[1]
				strand = info [2]
				start = info[3]
				end = info[4]
				dn = info[5]
				EST_coverage = info[6]
				cDNA_coverage = info[7]
				EST = info[8]			
				cDNA = info[9]						
		except KeyError:
			EST_coverage = 0
			cDNA_coverage = 0
			EST = 'NO'			
			cDNA = 'NO'
		
		
		try:
			info = bodymap_dict[DR_intron_ID]

			if bodymap_dict.has_key(DR_intron_ID)== True:
				intron = info[0]
				chr = info[1]
				strand = info [2]
				start = info[3]
				end = info[4]
				dn = info[5]		
				bodymap_coverage = info[6]
		except KeyError:
			bodymap_coverage = 0
			
		
		try:
			info = GM17828_dict[DR_intron_ID]
			
			if GM17828_dict.has_key(DR_intron_ID)== True:
				intron = info[0]
				chr = info[1]
				strand = info [2]
				start = info[3]
				end = info[4]
				dn = info[5]		
				GM17828_coverage = info[6]
		except KeyError:
			GM17828_coverage = 0		
		
	
		try:
			info = genecode_dict[DR_intron_ID]
			
			if genecode_dict.has_key(DR_intron_ID)== True:
				intron = info[0]
				chr = info[1]
				strand = info [2]
				start = info[3]
				end = info[4]
				dn = info[5]		
				genecode_transcript = info[6]
		except KeyError:
			genecode_transcript = 'NO'

	


#		if bodymap_coverage != 0: #and GM17828_coverage != 0:
		print intron, chr, strand, start, end, dn, bodymap_coverage, GM17828_coverage, genecode_transcript, EST_coverage, cDNA_coverage, EST, cDNA 




if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
