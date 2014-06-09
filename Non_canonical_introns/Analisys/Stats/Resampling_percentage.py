import sys
import csv
import pysam
import random
from collections import defaultdict

def main(final_table):
	
	reader1 = csv.reader(open(final_table), delimiter = ' ')
	coverage_list = []

	EST_cDNA = {}
	
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
		
		hg19_cDNA_EST = hg19_cDNA_coverage + hg19_EST_coverage
		mm9_cDNA_EST = mm9_cDNA_coverage + mm9_EST_coverage
		
		EST_cDNA[intron] = (hg19_cDNA_EST, mm9_cDNA_EST)  
		
		coverage = bodymap_coverage + gm12878_coverage
		
		for i in range(coverage):
			coverage_list.append((intron, dn, genecode_coverage))
		
		
	random.shuffle(coverage_list)  
	
	N_sample = int((len(coverage_list))/100)
	




	
	TOTAL = defaultdict(int)  #set([])
	canonical_known = defaultdict(int) #set([])
	canonical_unknown =  defaultdict(int) #set([])
	non_canonical_known =  defaultdict(int) #set([])
	non_canonical_unknown = defaultdict(int)  #set([])
	

	print "TOTAL", "canonical_known", "canonical_unknown", "non_canonical_known", "non_canonical_unknown" 
	
	
	for n in range(100):
		TOTAL_count = 0
		canonical_known_count = 0
		canonical_unknown_count =  0
		non_canonical_known_count =  0
		non_canonical_unknown_count = 0		
		
		sample_start = n * N_sample
		sample_end = (n+1) * N_sample
		if n + 1 == 100:
			for i in coverage_list[sample_start:]:
				intron = i[0]
				dn = i[1]
				genecode_coverage = i[2] 
				TOTAL[intron] += 1
				if dn == "GTAG" or dn == "GCAG" or dn == "ATAC":
					if genecode_coverage != 0:
						canonical_known[intron] += 1
					
					else:
						canonical_unknown[intron] += 1
				
				else:
					if genecode_coverage != 0:
						non_canonical_known[intron] += 1
					
					else:
						non_canonical_unknown[intron] += 1					
				
			
		else:
			for i in coverage_list[sample_start:sample_end]:
				intron = i[0]
				dn = i[1]
				genecode_coverage = i[2] 						
				TOTAL[intron] += 1
				if dn == "GTAG" or dn == "GCAG" or dn == "ATAC":
					if genecode_coverage != 0:
						canonical_known[intron] += 1
					
					else:
						canonical_unknown[intron] += 1
				
				else:
					if genecode_coverage != 0:
						non_canonical_known[intron] += 1
					
					else:
						non_canonical_unknown[intron] += 1	
		

		for i in TOTAL.items():
			intron = i[0]
			C = i[1]
			hg19_EST_cDNA = EST_cDNA[intron][0]
			mm9_EST_cDNA = EST_cDNA[intron][1]
			
			if  hg19_EST_cDNA>=3 or mm9_EST_cDNA>=3:			
				if C >=3:
					TOTAL_count += 1
			else:
				if C>=6:
					TOTAL_count += 1
					
		for i in canonical_known.items():
			intron = i[0]
			C = i[1]
			hg19_EST_cDNA = EST_cDNA[intron][0]
			mm9_EST_cDNA = EST_cDNA[intron][1]
			
			if  hg19_EST_cDNA>=3 or mm9_EST_cDNA>=3:	
				if C >=3:
					canonical_known_count += 1
			else:
				if C>=6:
					canonical_known_count += 1					 

		for i in canonical_unknown.items():
			intron = i[0]
			C = i[1]
			hg19_EST_cDNA = EST_cDNA[intron][0]
			mm9_EST_cDNA = EST_cDNA[intron][1]
			
			if  hg19_EST_cDNA>=3 or mm9_EST_cDNA>=3:	
				if C >=3:
					canonical_unknown_count += 1
			else:
				if C >=6:
					canonical_unknown_count += 1				

		for i in non_canonical_known.items():
			intron = i[0]
			C = i[1]
			hg19_EST_cDNA = EST_cDNA[intron][0]
			mm9_EST_cDNA = EST_cDNA[intron][1]
			
			if  hg19_EST_cDNA>=3 or mm9_EST_cDNA>=3:	
				if C >=3:
					non_canonical_known_count += 1
			else:
				if C >=6:
					non_canonical_known_count += 1

		for i in non_canonical_unknown.items():
			intron = i[0]
			C = i[1]
			hg19_EST_cDNA = EST_cDNA[intron][0]
			mm9_EST_cDNA = EST_cDNA[intron][1]
			
			if  hg19_EST_cDNA>=3 or mm9_EST_cDNA>=3:	
				if C >=3:
					non_canonical_unknown_count += 1
			else:
				if C >=6:
					non_canonical_unknown_count += 1				
		
		print TOTAL_count, canonical_known_count, canonical_unknown_count, non_canonical_known_count, non_canonical_unknown_count 

#		print len(TOTAL), len(canonical_known), len(canonical_unknown), len(non_canonical_known), len(non_canonical_unknown) 
		
		
		


if __name__ == '__main__':
	main(sys.argv[1])
