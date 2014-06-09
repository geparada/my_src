import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#USAGE : python ~/my_s

SeqTable = []

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


def main(pre_final_table, non_canonical_blat_genome, non_canonical_blat_tags):
		
	Genome = dict(SeqTable)	
	intron_blat_list = set([])
	intron_direct_blat_finded_list = set([])
	intron_bad_list = set([])
	intron_blat_not_found_list = set([])	
	intron_blat_tags = set([])
	
	introns_blat = set([])
	introns_mapsplice = set([])
	
	intron_info_list = []
	
	csv.field_size_limit(1000000000)
	
	reader1 = csv.reader(open(pre_final_table), delimiter = '\t')
	reader2 = csv.reader(open(non_canonical_blat_genome), delimiter = '\t')
	reader3 = csv.reader(open(non_canonical_blat_tags), delimiter = '\t')
	
	for row in reader3:
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
				
		block_tag_up = int(tName.split('|')[2].split('_')[0])
		block_tag_down = int(tName.split('|')[2].split('_')[1])
		
		intron = tName.split('|')[0]
				
		if misMatches == 0 and tStart<=(block_tag_up-8) and tEnd>=(block_tag_up+8) and (qEnd - qStart) > (qSize/2):     #revisar
			intron_blat_tags.add(intron)
	
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
		
		for b, t1, t2 in zip(blockSizes, tStarts, tStarts[1:]):
			
			ichr = tName
			istart = t1 + b
			iend = t2
			intron = tName + ':' + str(istart) + strand + str(iend)
			dn = Genome[ichr][istart:(istart+2)] + Genome[ichr][(iend-2):iend]
			
			if strand == '-':
				dn = dn.reverse_complement()
			
			dn = str(dn).upper()	
			
			DR_introns = DR_counter(intron, ichr, strand, istart, iend, Genome)
				
			DR_intron_ID = sorted(DR_introns)[0]
			
			introns_blat.add(str(dn)+DR_intron_ID)
			
			intron_blat_list.add((qName, (DR_intron_ID, dn)))
			
		
	intron_blat_dict = dict(intron_blat_list)
	
	for row in reader1:
		ichr = row[0]
		istart = int(row[1])
		iend = int(row[2])
		coverage = row[4]
		istrand = row[5]
		intron = ichr + ':' + str(istart) + istrand + str(iend)
		dn = Genome[ichr][istart:(istart+2)] + Genome[ichr][(iend-2):iend]        #Se extraen denuevo los dn para filtar todos los intrones que quedaron como canonicos en hg19
		if istrand == "-":
			dn = dn.reverse_complement()
		dn = str(dn).upper()
		reads = row[3].split("|")[2].split(',')
		ilength = iend - istart 
		
		intron_info_list.append((intron, [coverage, ichr, istrand, istart, iend, ilength, dn, ','.join(reads)]))   
		
		DR_introns = DR_counter(intron, ichr, istrand, istart, iend, Genome)
		
		if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
			print intron, coverage, ichr, istrand, istart, iend, ilength, dn, ','.join(reads)
		
		else:
			DR_intron_ID = sorted(DR_introns)[0]
			introns_mapsplice.add(str(dn)+DR_intron_ID)		
			
			for read in reads:
				
				try:
					dn_blat = intron_blat_dict[read][1]
					DR_intron_ID_blat = intron_blat_dict[read][0]
#					print read, DR_intron_ID, dn, DR_intron_ID_blat, dn_blat
					
					if DR_intron_ID == DR_intron_ID_blat:
						intron_direct_blat_finded_list.add(intron)
					
				except KeyError:
					intron_blat_not_found_list.add(intron)
					
	intron_info_dict = dict(intron_info_list)
					
	to_search_in_blat_tags = intron_blat_not_found_list    #Los intrones que BLAT los encontro como canonicos se descartan

	introns_suported_by_tags = to_search_in_blat_tags & intron_blat_tags   #Corresponden a los intrones que fueron rescatados con los tags
		
	intron_final_list = intron_direct_blat_finded_list | introns_suported_by_tags

#	print len(intron_final_list	)
#	print len(intron_final_list - intron_blat_tags)
	
	for intron in intron_final_list:
		info = intron_info_dict[intron] 
		print intron, " ".join(map(str,info))
		
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2],sys.argv[3],sys.argv[4])
