import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

		#matches = int(row[0])
		#misMatches = int(row[1])
		#repMaches = int(row[2])
		#nCount = int(row[3])
		#qNumInsert = int(row[4])    
		#qBaseInsert = int(row[5])
		#tNumInsert = int(row[6])
		#tBaseInsert = int(row[7])
		#strand = row[8]
		#qName = row[9]
		#qSize = int(row[10])
		#qStart = int(row[11])
		#qEnd = int(row[12])         
		#tName = row[13]
		#tSize = int(row[14])
		#tStart = int(row[15])
		#tEnd = int(row[16])
		#blockCount = int(row[17])
		#blockSizes = map(int, row[18].strip(',').split(','))
		#qStarts = map(int, row[19].strip(',').split(','))
		#tStarts = map(int, row[20].strip(',').split(','))

SeqTable = []
		
def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

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
				dn = Genome[ichr][istart+DRU:(istart+2+DRU)] + Genome[iend+DRU][(iend-2+DRU):iend]
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
		
def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0		

def main(SJ_introns, psl):
	Genome = dict(SeqTable)
	reader1 = csv.reader(open(psl), delimiter = '\t')
	reader2 = csv.reader(open(psl), delimiter = '\t')
	reader3 = csv.reader(open(SJ_introns), delimiter = '\t')
	reader4 = csv.reader(open(SJ_introns), delimiter = '\t')	
	
	
	alignment_scores = defaultdict(list)
	alignment_scores_dif = defaultdict(list)
	alignment_filtered = []
	
	blat = []
	bad_non_canonical = []
	
	for row in reader1:
		matches = int(row[0])
		misMatches = int(row[1])
		repMaches = int(row[2])
		qNumInsert = int(row[4])    
		tNumInsert = int(row[6])
		qName = row[9]

		
		score = matches + repMaches - misMatches - qNumInsert -tNumInsert
		alignment_scores[qName].append(score)
		
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
		
		score = matches + repMaches - misMatches - qNumInsert -tNumInsert	
		max_score = max(alignment_scores[qName])
		max_score_dif = percent(score, max_score)
					
		#for t, t2, b, q in zip(tStarts, tStarts[1:], blockSizes, qStarts[1:]):
		#	tbsum = t+b
		#	intron = tName + ":" + str(tbsum) + strand + str(t2)
		#	print qName, intron, score, max_score_dif
		
		if max_score_dif >= 99:
			alignment_scores_dif[qName].append(row)
	
	for alignment in alignment_scores_dif.items():          #convierte el dict en list y filtra los read que mapearon mas de una vez
		if len(alignment[1])==1:
			row = alignment[1][0]
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

			blat.append((qName, [tName,tStart,tEnd]))
			
			 

			#print '\t'.join(row)
	dict_blat = dict(blat)

			
	for row in reader3:

		read = row[0]
		chr = row[1]
		istart = row[2]
		iend = row[3]
		strand = row[4]
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
					   

		if dn!='GTAG' and dn!='GCAG' and dn!='ATAC':
			try:
				DR_canonical = DR_Canonical_intron_finder(intron, chr, strand, int(istart), int(iend), Genome)
				
				if len(DR_canonical)>0:
					bad_non_canonical.append((read, [intron, dn, 'Hidden_DR_Canonical'] ))					
			
				blat_chr = dict_blat[read][0]
				blat_start = dict_blat[read][1]
				blat_end = dict_blat[read][2]
			
				if chr == blat_chr and ((start <= blat_start < end) or (blat_start <= start < blat_end)):
					pass
					
				else:
					bad_non_canonical.append((read, [intron, dn, 'Wrong'] ))						
					
			except KeyError:
				bad_non_canonical.append((read, [intron, dn, 'Multimapping'] ))
	
	dict_bad_non_canonical = dict(bad_non_canonical)
	
	for row in reader4:
		read = row[0]
		
		if dict_bad_non_canonical.has_key(read) == False:
			print '\t'.join(row)
			
		else:
			print '\t'.join(row), '\t' + 'BAD_ALIGNMENT' + '\t', '\t'.join(dict_bad_non_canonical[read])
			
			
# ******Inicio******

#La idea es que este programa lea el SJ.bed12, y si es que hay un intron no canonico, entonses vea el psl de BLAT
#Se crea un dicionario (dict_blat) de todos los intrones no canonicos que blat encontron en reads que alineaban solo una vez al genoma

#*****Clasificacion de errores****

#Si no se en cuentra el key en dict_blat ---> se clasifica como Multimapping (aun que tambien pudo haber sido porque blat no pudo alinear el read)
#Si los alieamientos de MapSplice y BLAT no se solapan ----> se clasifica como Wrong
#Si al correr el alineamiento por sus directos repetidos se encuentra un intron canonico -----> se clasifica como Hidden_DR_Canonical


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3])  
	
	#USAGE: python BLAT_filter.py genome.fa bed12 psl
