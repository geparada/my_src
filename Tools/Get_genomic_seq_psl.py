import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from operator import itemgetter, attrgetter


Genome = {}


def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[str(chrfa.id)] =  str(chrfa.seq)

	f.close()
	


def main(psl):
	
	reader = csv.reader(open(psl), delimiter = '\t')
	
	results = {}
	IDs = []
	
	for row in reader:
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
		
		seq = Genome[tName][tStart-10:tEnd+10]
		
		if strand=="-":
			seq = str(Seq(seq).reverse_complement())
		
		ID = qName.split("|")[2] + "|" + tName.split("|")[3]
		
		if (ID in IDs) == False:
			IDs.append(ID)
			
		score = matches + repMaches - misMatches - qNumInsert
		
		if results.has_key(ID)==False:
			results[ID] = seq, score
		
		elif results[ID][1] < score:
			results[ID] = seq, score
	
	
	result_list = sorted(results.items(), key=itemgetter(0))
	
		
	for i in result_list:
		if "hugZ" in i[0]:
			print ">" + i[0]
			print i[1][0]



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])
			
