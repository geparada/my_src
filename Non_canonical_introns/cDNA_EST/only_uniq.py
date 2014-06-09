import sys
import csv
from collections import defaultdict


def main(psl):
	
	reader = csv.reader(open(psl), delimiter = '\t')
	alignments = defaultdict(list)
	
	for row in reader:

		matches = int(row[1])
		misMatches = int(row[2])
		repMaches = int(row[3])
		nCount = int(row[4])
		qNumInsert = int(row[5])    
		qBaseInsert = int(row[6])
		tNumInsert = int(row[7])
		tBaseInsert = int(row[8])
		strand = row[9]
		qName = row[10]
		qSize = int(row[11])
		qStart = int(row[12])
		qEnd = int(row[13])         
		tName = row[14]
		tSize = int(row[15])
		tStart = int(row[16])
		tEnd = int(row[17])
		blockCount = int(row[18])
		blockSizes = map(int, row[19].strip(',').split(','))
		qStarts = map(int, row[20].strip(',').split(','))
		tStarts = map(int, row[21].strip(',').split(','))
		
		alignments[qName].append(row[1:])
	
	for row in alignments.items():
		qName = row[0]
		alignments = row[1]
		if len(alignments) == 1:
			print "\t".join(alignments[0])
		
		
		
		
		
		
		
if __name__ == '__main__':
	main(sys.argv[1])
