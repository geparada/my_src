import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict

Genome = {}
Transcriptome = {}
out1 = open("tags", "w")

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()
	
def Transcriptabulator(fasta):
	
	print >> sys.stderr, "Cargando cDNAs/ESTs en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Transcriptome[chrfa.id] = chrfa.seq

	print >> sys.stderr, "OK"

	f.close() 

def DR_Canonical_intron_finder(intron, ichr, strand, istart, iend, Genome):
	
	'''Esta funcion busca intrones canonicos (GTAG, ATAC, GCAG) al correr el intron por sus directos repetidos''' 
	
	introns_finded_DR = []
			
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
				dn = Genome[ichr][(istart-DRU):(istart+2-DRU)] + Genome[ichr][(iend-2-DRU):iend-DRU]
				dn = str(dn).upper()
				if dn == "GTAG" or dn == "GCAG" or dn == "ATAC": 								
					introns_finded_DR.append(dn)
			elif strand == "-":
				dn = Genome[ichr][istart+DRU:(istart+2+DRU)] + Genome[ichr][(iend-2+DRU):iend+DRU]
				dn = str(dn.reverse_complement()).upper()
				if dn == "GTAG" or dn == "GCAG" or dn == "ATAC": 								
					introns_finded_DR.append(dn)
			if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
				break
	except IndexError:
		pass 
	try:
		while SJ5D[DRD]==SJ3D[DRD]:
			DRD += 1
			if strand == "+":
				dn = Genome[ichr][(istart+DRD):(istart+2+DRD)] + Genome[ichr][(iend-2+DRD):iend+DRD]
				dn = str(dn).upper()				
				if dn == "GTAG" or dn == "GCAG" or dn == "ATAC": 								
					introns_finded_DR.append(dn)
			elif strand == "-":
				dn = Genome[ichr][istart-DRD:(istart+2-DRD)] + Genome[ichr][(iend-2-DRD):iend-DRD]
				dn = str(dn.reverse_complement()).upper()
				if dn == "GTAG" or dn == "GCAG" or dn == "ATAC": 								
					introns_finded_DR.append(dn)				
			if SJ5D[DRD]!=SJ3D[DRD]:
				break
	except IndexError:
		pass
	
	return introns_finded_DR

	
def main(psl):
	
	reader2 = csv.reader(open(psl), delimiter = '\t')
	
	for row in reader2:
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

		dn_count = defaultdict(int)
				
		
		
		for b, t1, t2 in zip(blockSizes, tStarts, tStarts[1:]):
			
			chr = tName
			istart = t1 + b
			iend = t2
			ilength = iend - istart
			intron = tName + ':' + str(istart) + strand + str(iend)
			dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
			
			if strand == '-':
				dn = dn.reverse_complement()
				
			if ilength >=40:				
				dn = str(dn).upper()
				dn_count[dn] += 1
	
			
		orient = 1
		if strand == '-':
			orient = -1		
			
		if ((dn_count['GTAG'] + dn_count['GCAG'] + dn_count['ATAC'])  - (dn_count['CTAC'] + dn_count['CTGC'] + dn_count['GTAT'] )) < 0:
			orient = -1 * orient
		
		fixed_strand = ''
		
		if orient == 1:
			fixed_strand = "+"
		elif orient == -1:
			fixed_strand = "-"
			
		
		for b, t1, t2, q1, q2 in zip(blockSizes, tStarts, tStarts[1:], qStarts, qStarts[1:]):
			
			chr = tName
			istart = t1 + b
			iend = t2
			ilength = iend - istart
			intron = tName + ':' + str(istart) + fixed_strand + str(iend)
			dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
			
			Qtag = Transcriptome[qName][(q1 + b - 15):(q1 + b)]  + Transcriptome[qName][q2: (q2 + 15)]			
			Gtag = Genome[chr][(istart-15):istart] + Genome[chr][iend:(iend+15)]
			
			if strand == '-':
				Qtag = Transcriptome[qName].reverse_complement()[(q1 + b - 15):(q1 + b)]  + Transcriptome[qName].reverse_complement()[q2: (q2 + 15)]	
			
			Qtag = str(Qtag).upper()
			Gtag = str(Gtag).upper()	
			

			out1.write(">" + qName + "|" + intron + "\n")
			out1.write(Qtag + "\n")
			
			if fixed_strand == '-':
				dn = dn.reverse_complement()
				
			dn = str(dn).upper()

			if ilength > 6 and Qtag == Gtag:   #Los gaps de 6 o menos son conciderados como deleciones por MapSplice
			
				if dn == "GTAG" or dn == "GCAG" or dn == "ATAC":
					print qName, chr, istart, iend, fixed_strand, ilength, intron, dn

				
				else:
					DR_canonical = DR_Canonical_intron_finder(intron, chr, fixed_strand, int(istart), int(iend), Genome)
					if len(DR_canonical) == 0:
						print qName, chr, istart, iend, fixed_strand, ilength, intron, dn
			
				
			
			
				



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	Transcriptabulator(sys.argv[2])
	main(sys.argv[3])
