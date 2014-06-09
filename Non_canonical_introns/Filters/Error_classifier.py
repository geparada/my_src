import sys
import csv
from collections import defaultdict

Genome = {}
multimaping_reads = set([])
read_ID_intron_DR_BLAT = set([])
read_ID_intron_DR_MapSplice = set([])
BLAT_info = {}

def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[str(chrfa.id)] = chrfa.seq

	f.close()


def dn_introns_DR_extractor(intron, ichr, strand, istart, iend, Genome):
	
	'''Esta funcion entrega:
	 -Los dinucleotidos de todos los posibles intrones 
	 -Un ID del intron basado en sus DR, necesario para comparar los intrones MapSplice/BLAT''' 
	

	istart_DR = [istart]
	iend_DR = [iend]
	
	dn_introns_DR = []
	introns_DR = []
			
	L = 100         #Solo permite que se corra L pares de bases para buscar DR 		

	#Extrayendo regiones exonicas colindantes
				
	SJ5U = Genome[ichr][istart-L : istart].lower()
	SJ5D = Genome[ichr][istart : istart+L].lower()
	SJ3U = Genome[ichr][iend-L : iend].lower()
	SJ3D = Genome[ichr][iend : iend+L].lower()
	
	dn = Genome[ichr][istart:istart+2] + Genome[ichr][iend-2:iend]
	
	if strand == "-":	
		dn = dn.reverse_complement()

	
	dn = str(dn).upper()
	dn_introns_DR.append(dn)
	
		
	DRU = 0
	DRD = 0
								
	#Contando directos repetidos y generando intrones no consenso alternativos
		
	try:
		while SJ5U[L-1-DRU]==SJ3U[L-1-DRU]:
			DRU += 1
			if strand == "+":
				
				dn = Genome[ichr][istart-DRU:istart+2-DRU] + Genome[ichr][iend-2-DRU:iend-DRU]
				dn = str(dn).upper()
				dn_introns_DR.append(dn)
				
				introns_DR.append(ichr + ":" + str(istart-DRU) + strand + str(iend-DRU))
				
				istart_DR.append(istart-DRU)
				iend_DR.append(iend-DRU)
				
			elif strand == "-":
				
				dn = Genome[ichr][istart+DRU:istart+2+DRU] + Genome[ichr][iend-2+DRU:iend+DRU]
				dn = str(dn.reverse_complement()).upper()
				dn_introns_DR.append(dn)
				
				introns_DR.append(ichr + ":" + str(istart-DRU) + strand + str(iend-DRU))				

				istart_DR.append(istart+DRU)
				iend_DR.append(iend+DRU)								

			if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
				break
	except IndexError:
		pass
		 
	try:
		while SJ5D[DRD]==SJ3D[DRD]:
			DRD += 1
			if strand == "+":
				
				dn = Genome[ichr][istart+DRD:istart+2+DRD] + Genome[ichr][iend-2+DRD:iend+DRD]
				dn = str(dn).upper()
				dn_introns_DR.append(dn)				
				
				introns_DR.append(ichr + ":" + str(istart-DRU) + strand + str(iend-DRU))
				
				istart_DR.append(istart+DRD)
				iend_DR.append(iend+DRD)
				
			elif strand == "-":
				
				dn = Genome[ichr][istart-DRD:istart+2-DRD] + Genome[ichr][iend-2-DRD:iend-DRD]
				dn = str(dn.reverse_complement()).upper()
				dn_introns_DR.append(dn)
				
				introns_DR.append(ichr + ":" + str(istart-DRU) + strand + str(iend-DRU))				
				
				istart_DR.append(istart-DRD)
				iend_DR.append(iend-DRD)
								
			if SJ5D[DRD]!=SJ3D[DRD]:
				break
	except IndexError:
		pass
	
	return dn_introns_DR, ",".join(sorted(dn_introns_DR))

def BLAT_ERROR ():
	""" Analiza el archivo pls para buscar multimaping y genera dicionario con lugares donde alinearon los reads """
	
	reader1 = csv.reader(open(psl), delimiter = '\t')
	reader2 = csv.reader(open(psl), delimiter = '\t')
	
	score_aligments_read = defaultdict(list)
	high_score_aligments = defaultdict(list)
	
	for row in reader1:
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
		
		score = matches + repMaches - misMatches - qNumInsert
		score_aligments_read[qName].append(score)

	
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
		
		score = matches + repMaches - misMatches - qNumInsert
		max_score = max(score_aligments_read[qName])
		
		if max_score - score < max_score/100:
			high_score_aligments[read].append(row)
		
		for b, t1, t2 in zip(blockSizes, tStarts, tStarts[1:]):
			
			ichr = tName
			istart = t1 + b
			iend = t2
			intron = tName + ':' + str(istart) + strand + str(iend)
			dn = Genome[ichr][istart:(istart+2)] + Genome[ichr][(iend-2):iend]
			
			if strand == '-':
				dn = dn.reverse_complement()
				
			dn_introns_DR, ID_introns_DR = dn_introns_DR_extractor(intron, ichr, pair_strand, istart, iend, Genome)
			
			read_ID_intron_DR_BLAT.add((qName, ID_introns_DR))
		
	
	for a in high_score_aligments.items():
		if len(a[1]>1):
			multimaping_reads.add(a[0])
		else:
			BLAT_info[qName] = tName, tStart, tEnd 
			 


def main ():
	""" Function doc """
	
	introns_MapSplice_BLAT = set([])  #Contiene todos los ID_DRs de los intrones que se encontraron por blat y por MapSplice
	
	reader1 = csv.reader(open(SJ_introns), delimiter = ' ')
	reader2 = csv.reader(open(SJ_introns), delimiter = ' ')
	
	for row in reader1:
		read = row[0]
		chr = row[1]
		istart = int(row[2])
		iend = int(row[3])
		pair_strand = row[4]
		ilen = row[5]
		intron = row[6]
		dn = row[7]
		start = row[8]
		cigar = row[9]
		e5s = row[10]
		e5e = row[11]
		e3s = row[12]
		e3e = row[13]
		seq = row[14]
		
		dn_introns_DR, ID_introns_DR = dn_introns_DR_extractor(intron, ichr, pair_strand, istart, iend, Genome)
				
		read_ID_intron_DR_MapSplice.add((read, ID_introns_DR))
	
		for i in read_ID_intron_DR_MapSplice & read_ID_intron_DR_MapSplice:
			introns_MapSplice_BLAT.add(i[1])


	for row in reader2:
		read = row[0]
		chr = row[1]
		istart = int(row[2])
		iend = int(row[3])
		pair_strand = row[4]
		ilen = row[5]
		intron = row[6]
		dn = row[7]
		start = int(row[8])
		cigar = row[9]
		e5s = row[10]
		e5e = row[11]
		e3s = row[12]
		e3e = row[13]
		seq = row[14]
		
		dn_introns_DR, ID_introns_DR = dn_introns_DR_extractor(intron, ichr, pair_strand, istart, iend, Genome)
		
		if read in multimaping_reads:
			print " ".join(row), BAD_ALIGNMENT, BLAT_multimaping
		
		else:
			BLAT_chr, BLAT_start, BLAT_end = BLAT_info[read]
		
			if "GTAG" in dn_introns_DR or "GCAG" in dn_introns_DR or "ATAC" in dn_introns_DR:
				print " ".join(row), BAD_ALIGNMENT, BLAT_canonical
			
			elif chr != BLAT_chr or start > BLAT_end or BLAT_end > (start + len(seq)):
				print " ".join(row), BAD_ALIGNMENT, DIFERENT_LOCUS
			
			elif intron in introns_MapSplice_BLAT:
				print " ".join(row)	
	
	


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	BLAT_ERROR(sys.argv[3])
	main(sys.argv[2])


#Se concidera multimaping cuando hay dos alinemientos con score maximo
#Debido al el softcpliping, basta que uno de los reads tenga el intron no canonico 	
