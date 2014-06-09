import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict

Genome = {}
out1 = open("psl_uniq", "w")

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
	
	print >> sys.stderr, "OK"

	f.close() 

def main(psl):
	reader1 = csv.reader(open(psl), delimiter = '\t')
	
	scores = defaultdict(list)
	scores_info = {}
	blat = {}	
	multimapping = {}
	
		
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
		qName = row[9].split("|")[0]
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
		
		score = matches + repMaches - repMaches - qNumInsert
		scores[qName].append(score)
		scores_info[qName+str(score)] = row  
	
	for row in scores.items():
		ID = row[0]
		max_score = max(row[1])
		max_score_psl = scores_info[ID+str(max_score)]
		blat[ID] = max_score_psl
		 
		high_relative_scores = []
		for s in row[1]:
			if s>= 99:
				high_relative_scores.append(percent(s, max_score))
		
		if len(high_relative_scores) > 1:
			multimapping[ID] =  len(high_relative_scores)
	
	for row in blat.items():
		matches = int(row[1][0])
		misMatches = int(row[1][1])
		repMaches = int(row[1][2])
		nCount = int(row[1][3])
		qNumInsert = int(row[1][4])    
		qBaseInsert = int(row[1][5])
		tNumInsert = int(row[1][6])
		tBaseInsert = int(row[1][7])
		strand = row[1][8]
		qName = row[1][9].split("|")[0]
		qSize = int(row[1][10])
		qStart = int(row[1][11])
		qEnd = int(row[1][12])         
		tName = row[1][13]
		tSize = int(row[1][14])
		tStart = int(row[1][15])
		tEnd = int(row[1][16])
		blockCount = int(row[1][17])
		blockSizes = map(int, row[1][18].strip(',').split(','))
		qStarts = map(int, row[1][19].strip(',').split(','))
		tStarts = map(int, row[1][20].strip(',').split(','))
		
		if multimapping.has_key(qName)==False:

			out1.write('\t'.join(row[1][:9]) + '\t' + qName  + '\t' + '\t'.join(row[1][10:]) + "\n")
			transcript_seq = ""
			
			for b, t1 in zip(blockSizes, tStarts):
					
				chr = tName
				exon = Genome[chr][t1:t1+b]
				transcript_seq += exon

			if strand == '-':
				transcript_seq = transcript_seq.reverse_complement()
			transcript_seq = str(transcript_seq).upper()
						
			print ">" + qName 
			print transcript_seq
			
						
			
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])
