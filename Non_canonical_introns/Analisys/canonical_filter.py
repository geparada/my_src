import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def main(Final_table):
	
	reader1 = csv.reader(open(Final_table), delimiter = ' ')
	
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

		introns = set([])

		if dn =="GTAG" or dn == "GCAG" or dn == "ATAC":
		
			info = " ".join(row)
			
			if bodymap_coverage>=3 and (hg19_cDNA_coverage + hg19_EST_coverage)>=3: 
				introns.add(info)
				
			
			if gm12878_coverage>=3 and (hg19_cDNA_coverage + hg19_EST_coverage)>=3: 
				introns.add(info)
			
			if bodymap_coverage>=3 and gm12878_coverage>=3: 		
				introns.add(info)
			
			if (bodymap_coverage>=3 and (mm9_cDNA_coverage + mm9_EST_coverage)>=3)  or (gm12878_coverage>=3 and (mm9_cDNA_coverage + mm9_EST_coverage)>=3):
				introns.add(info)
		
		
		for i in introns:
			print i, "*", "*", "*"

#python ~/my_src/Analisys/SNPs_filter.py ~/db/genome/hg19.fa ~/db/Variation/snp135.txt non_canonical.final_table.tags.seq.DR

if __name__ == '__main__':
	main(sys.argv[1])
