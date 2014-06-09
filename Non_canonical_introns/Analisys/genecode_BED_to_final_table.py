import sys
import csv
from collections import defaultdict


def main(BED):
	reader1 = csv.reader(open(BED), delimiter = ' ')
	reader2 = csv.reader(open(BED), delimiter = ' ')
	
	genecode_introns = defaultdict(list)
	
	for row in reader1:
		transcript = row[0]
		chr = row[1]
		istart = row[2]
		iend = row[3]
		strand = row[4]
		ilenght = row[5]
		intron = row[6]
		dn = row[7]
		
		genecode_introns[intron].append(transcript)
		
	for row in reader2:
		chr = row[1]
		istart = row[2]
		iend = row[3]
		strand = row[4]
		ilenght = row[5]
		intron = row[6]
		dn = row[7]
		
		transcripts = genecode_introns[intron]
		
		print intron, len(transcripts), chr, strand, istart, iend, ilenght, dn, ",".join(transcripts)	
		
	



if __name__ == '__main__':
	main(sys.argv[1])
