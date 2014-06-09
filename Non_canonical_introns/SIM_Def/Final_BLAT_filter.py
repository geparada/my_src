import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable=[]


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

	f.close()  

	print >> sys.stderr, "OK"

def main(introns_final, non_can_introns, non_cannonical_psl):
	Genome = dict(SeqTable)
	
	BLAT_psl = defaultdict(int)
	BLAT_reads_introns = defaultdict(list)
	
	Final_introns_BLAT = defaultdict(set)
	
	non_can_intron_read = defaultdict(list)
	non_can_list = []
	
	reader1 = csv.reader(open(introns_final), delimiter = ' ')
	reader2 = csv.reader(open(non_cannonical_psl), delimiter =  '\t')
	reader3 = csv.reader(open(non_can_introns), delimiter = ' ')
	reader4 = csv.reader(open(non_cannonical_psl), delimiter =  '\t')
	reader5 = csv.reader(open(introns_final), delimiter = ' ')
	
	print >> sys.stderr, "Extrayendo nombres de reads que tienen intrones no canonicos ...",
	
	for row in reader3:
		read = row[0]
		intron = row[6]
		non_can_intron_read[intron].append(read)
		
	print >> sys.stderr, "OK"
	
	print >> sys.stderr, "Extrayendo intrones no canonicos de reads que solo alinearon una vez con BLAT ...",
	
	for row in reader2:

		qName = row[9]
		BLAT_psl[qName] += 1
			
	for row in reader4:
		matches = int(row[0])
		misMatches = int(row[1])
		repMaches = int(row[2])
		nCount = int(row[3])
		qNumInsert = int(row[4])    
		qBaseInsert = int(row[5])
		tNumInsert = int(row[6])
		tBaseInsert = int(row[7])
		strand = row[8]
		qName = row[9]
		qSize = int(row[10])
		qStart = int(row[11])
		qEnd = int(row[12])         
		tName = row[13]
		tSize = int(row[14])
		tStart = int(row[15])
		tEnd = int(row[16])
		blockCount = int(row[17])
		blockSizes = map(int, row[18].strip(',').split(','))
		qStarts = map(int, row[19].strip(',').split(','))
		tStarts = map(int, row[20].strip(',').split(','))
		
		if BLAT_psl[qName]==1:             #Solo se toman en cuenta los intrones de los reads que alinearon una vez
			for b, t1, t2 in zip(blockSizes, tStarts, tStarts[1:]):
				intron = tName + ':' + str(t1 + b) + strand + str(t2)
				BLAT_reads_introns[qName].append(intron)
				
	print >> sys.stderr, "OK"
	
	print >> sys.stderr, "Buscando concidencias de intrones no canonicos con BLAT ...",
	
	for row in reader1:
		intron = row[0]
		strand = row[3]
		dn = row[7]
		ichr = row[2]
		istart = int(row[4])
		iend = int(row[5])
		
		if dn!='GTAG' and dn!='GCAG' and dn!='ATAC':
			reads = non_can_intron_read[intron]
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
						
			intron_reads = non_can_intron_read[intron]
			for read in intron_reads:
				for blat_intron in BLAT_reads_introns[read]:
					if blat_intron in introns_finded_DR:
						Final_introns_BLAT[introns_finded_DR[0]].add("OK")    #El primer intron de la lista introns_finded_DR es el que corresponde al encontrado por MapSplice/GSNAP				
	
	print >> sys.stderr, "OK"
	
	print >> sys.stderr, "Aplicando filtro ...",
	
	for row in reader5:
		intron = row[0]
		dn = row[7]
		
		if dn=='GTAG' or dn=='GCAG' or dn=='ATAC':
			print " ".join(row)
			
		elif Final_introns_BLAT.has_key(intron)==True:
			print " ".join(row)

		
	print >> sys.stderr, "OK"			
			
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2],sys.argv[3],sys.argv[4]) 	
	
			
			
		
	
