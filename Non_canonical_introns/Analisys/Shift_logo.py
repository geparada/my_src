import sys
import csv
from collections import defaultdict
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
	
	return DRU, DRD

def main (FINAL_TABLE):
	""" Extrae secuencias par logo de los shift """
	
	reader1 = csv.reader(open(FINAL_TABLE), delimiter="\t")
	
	for row in reader1:
		
		intron = row[1]
		strand = row[3]
		shifts = row[26].split(",")
		for i in shifts:
				alt_intron, dn, fold, donor, aceptor = i.split("|")
				donor = int(donor)
				aceptor = int(aceptor)
				
				ichr = alt_intron.split(":")[0]
				istart = int(alt_intron.split(":")[1].replace("+","-").split("-")[0])
				iend = int(alt_intron.split(":")[1].replace("+","-").split("-")[1])

				
				DRU, DRD = DR_counter(alt_intron, ichr, strand, istart, iend, Genome)
		
				print intron, i, DRU, DRD
				
#				for n in range(DRU):
					
	
	



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2])
