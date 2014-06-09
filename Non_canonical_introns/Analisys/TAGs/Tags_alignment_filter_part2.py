import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict

Genome = {}

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

def main(tags_filter, reads_genome_blat, reads_repbase, reads_dust ):
	reader1 = csv.reader(open(tags_filter), delimiter = ' ')	
	reader2 = csv.reader(open(reads_genome_blat), delimiter = '\t')
	reader3 = csv.reader(open(reads_repbase), delimiter = '\t')
	reader4 = csv.reader(open(reads_dust), delimiter = ' ')


	blat = {}	
	repbase = {}
	dust = {}
	multimapping = {}
	
	
	for row in reader3:
		read = row[9]
		repeat = row[13]
		repbase[read] = repeat
	
	for row in reader4:
		read = row[0].strip(">")
		dust[read] = 0
	
	scores = defaultdict(list)
	scores_info = {}
		
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
		
		score = matches + repMaches - repMaches - qNumInsert
		scores[qName].append(score)
		scores_info[qName+str(score)] = row  
	
	for row in scores.items():
		read = row[0]
		max_score = max(row[1])
		max_score_psl = scores_info[read+str(max_score)]
		blat[read] = max_score_psl
		 
		high_relative_scores = []
		for s in row[1]:
			if s>= 99:
				high_relative_scores.append(percent(s, max_score))
		
		if len(high_relative_scores) > 1:
			multimapping[read] =  len(high_relative_scores)		
	
	for row in reader1:
		read = row[0]
		seq = row[1]
		qual = row[2]
		intron = row[3]
		up_anchor = int(row[4])
		down_anchor = int(row[5])
		
		try:
		
			matches = int(blat[read][0])
			misMatches = int(blat[read][1])
			repMaches = int(blat[read][2])
			nCount = int(blat[read][3])
			qNumInsert = int(blat[read][4])    
			qBaseInsert = int(blat[read][5])
			tNumInsert = int(blat[read][6])
			tBaseInsert = int(blat[read][7])
			strand = blat[read][8]
			qName = blat[read][9]
			qSize = int(blat[read][10])
			qStart = int(blat[read][11])
			qEnd = int(blat[read][12])         
			tName = blat[read][13]
			tSize = int(blat[read][14])
			tStart = int(blat[read][15])
			tEnd = int(blat[read][16])
			blockCount = int(blat[read][17])
			blockSizes = map(int, blat[read][18].strip(',').split(','))
			qStarts = map(int, blat[read][19].strip(',').split(','))
			tStarts = map(int, blat[read][20].strip(',').split(','))
			
			if qSize - (qEnd - qStart) > min(up_anchor, down_anchor):
				if repbase.has_key(read)==False and dust.has_key(read)==False and multimapping.has_key(read)==False:
					print " ".join(row)
			
			else:
				non_canonical_introns = []
				
				for b, t1, t2 in zip(blockSizes, tStarts, tStarts[1:]):
				
					ichr = tName
					istart = t1 + b
					iend = t2
					intron = tName + ':' + str(istart) + strand + str(iend)
					dn = Genome[ichr][istart:(istart+2)] + Genome[ichr][(iend-2):iend]
					
					if strand == '-':
						dn = dn.reverse_complement()
					
					dn = str(dn).upper()
					
					if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
						non_canonical_introns.append(dn)
				
				if len(non_canonical_introns)!=0 and repbase.has_key(read)==False and dust.has_key(read)==False and multimapping.has_key(read)==False:
					print " ".join(row)
						
		
		except KeyError:
			
			if repbase.has_key(read)==False and dust.has_key(read)==False and multimapping.has_key(read)==False:
				print " ".join(row)
		

		
	 
	
	
		
	
	
	

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
