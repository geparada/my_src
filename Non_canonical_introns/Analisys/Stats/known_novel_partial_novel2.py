import sys
import csv
from collections import defaultdict
from decimal import *


def main(genecode_introns, final_table):

	reader1 = csv.reader(open(genecode_introns), delimiter = ' ')		
	reader2 = csv.reader(open(final_table), delimiter = ' ')
	

	genecode_introns = set([])
	
	for row in reader1:
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

		alt_introns = row[21].split(",")
		
		if genecode_coverage != 0:
			genecode_introns.add(intron)



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

#		intron_retention_exon = row[19].split(",")
#		skipped_exons_names = row[20].split(",")
		
		alt_introns = row[21].split(",")
		
#		alt_no_skipper_introns = row[22].split(",")
#		alt_skipper_introns = row[23].split(",")
#		alt_exon_variant_introns = row[24].split(",")
#		shift = row[25].split(",")
#		non_canonical_filter = row[26].split(",")		
		
		
		if genecode_coverage!=0:
			known_TOTAL += 1
			if dn == "GTAG" or dn == "GCAG" or dn == "ATAC":
				known_canonical += 1				
			else:
				known_non_canonical += 1			
						
		
		k5 = chr+strand+start
		k3 = chr+strand+end
		
		elif k5 in genecode_5 or k3 in genecode_3:
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
				print row 
				novel_non_canonical += 1					
		
	print known_TOTAL, partial_novel_TOTAL, novel_TOTAL
	print known_canonical, partial_novel_canonical, novel_canonical		
	print known_non_canonical, partial_novel_non_canonical, novel_non_canonical			
			 		


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])


