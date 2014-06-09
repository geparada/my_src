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
		Genome[chrfa.id] = chrfa.seq

	f.close() 
	
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
		
		for b, t1, t2 in zip(blockSizes, tStarts, tStarts[1:]):
			
			ichr = tName
			istart = t1 + b
			iend = t2
			ilength = iend - istart
			intron = tName + ':' + str(istart) + strand + str(iend)
			dn = Genome[ichr][istart:(istart+2)] + Genome[ichr][(iend-2):iend]
			
			if strand == '-':
				dn = dn.reverse_complement()
				
			print qName, chr, istart, iend, strand, ilength, intron, dn


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])

