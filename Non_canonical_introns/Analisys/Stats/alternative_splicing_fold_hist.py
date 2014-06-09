import sys
import csv
from decimal import *
import math
getcontext().prec = 20
from collections import defaultdict


def log_fold(fold):
	if fold == 0:
		return 0
	else:
		return math.log(fold**-1,2)			

def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')
	fold_histogram = defaultdict(int)

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
		intron_retention_exon = row[19].split(",")
		skipped_exons_names = row[20].split(",")
		alt_introns = row[21].split(",")
		alt_no_skipper_introns = row[22].split(",")
		alt_skipper_introns = row[23].split(",")
		alt_exon_variant_introns = row[24].split(",")
		shift = row[25].split(",")
		non_canonical_filter = row[26].split(",")
		
		folds = []
		
		if intron_retention_exon!=['NO']:
			for i in intron_retention_exon:
				folds.append( Decimal(i.split("|")[1]))

		if skipped_exons_names!=['NO']:
			for i in skipped_exons_names:
				folds.append( Decimal(i.split("|")[1]))

		if alt_introns!=['NO']:
			for i in alt_introns:
				folds.append( Decimal(i.split("|")[2]))

		if shift!=['NO']:
			for i in shift:
				folds.append( Decimal(i.split("|")[2]))
				
		max_fold = 0
		
		if len(folds)!=0:
			max_fold = max(folds)

		if max_fold!=0:
			fold_histogram[int(log_fold(max_fold))] += 1

	for n in range(-14,17):
		print str(n - 1) +"\t" + str(fold_histogram[n])
	 

			
			
if __name__ == '__main__':
	main(sys.argv[1])
