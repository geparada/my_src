import sys
import csv
from collections import defaultdict
from decimal import *


def main(genecode_introns, final_table):

	reader1 = csv.reader(open(genecode_introns), delimiter = ' ')		
	reader2 = csv.reader(open(final_table), delimiter = ' ')
	

	genecode_5 = set([])
	genecode_3 = set([]) 
	
	
	for row in reader1:
		ID = row[0]
		chr = row[1]
		start = row[2]
		end = row[3]
		strand = row[4]
		lenght = row[5]
		intron = row[6]
		dn = row[7]
		
		genecode_5.add(chr+strand+start)
		genecode_3.add(chr+strand+end)


	known_TOTAL = 0
	partial_novel_TOTAL = 0
	novel_TOTAL = 0
	
	known_canonical = 0
	partial_novel_canonical = 0
	novel_canonical = 0

	known_non_canonical = 0
	partial_novel_non_canonical = 0
	novel_non_canonical = 0

	
	for row in reader2:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = row[8]
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		bodymap_seq = row[16]
		gm12878_seq = row[17]
		DR = row[18]

		
		if genecode_coverage!=0:
			known_TOTAL += 1
			if dn == "GTAG" or dn == "GCAG" or dn == "ATAC":
				known_canonical += 1				
			else:
				known_non_canonical += 1			

		
		elif chr+strand+str(istart) in genecode_5 or chr+strand+str(iend) in genecode_3:
			partial_novel_TOTAL += 1
			if dn == "GTAG" or dn == "GCAG" or dn == "ATAC":
				partial_novel_canonical += 1				
			else:
				partial_novel_non_canonical += 1				 		

		else:
			novel_TOTAL += 1
			if dn == "GTAG" or dn == "GCAG" or dn == "ATAC":
				novel_canonical += 1				
			else:
#				print intron, dn, dn_type_score, row[21] 
				novel_non_canonical += 1					
		
	print known_TOTAL, partial_novel_TOTAL, novel_TOTAL
	print known_canonical, partial_novel_canonical, novel_canonical		
	print known_non_canonical, partial_novel_non_canonical, novel_non_canonical			
			 		


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])


