import sys
import csv
from collections import defaultdict

#Se debe tomar el la lista que tambien tiene los canonicos!

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	



def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')
	
	

	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)
	
	U2_GTAG = []
	U12_ATAC = []
	U12_GTAG = []
	indef = []

	dn_total = defaultdict(int)
	dn_U2 = defaultdict(int)
	dn_U12 = defaultdict(int)	


	for row in reader1:
		row = row[1:]

		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = float(row[8])
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
		
		dn_total[dn] += 1
		

		
		if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
			U2_GTAG.append(row)
		
		elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
			U12_ATAC.append(row)		
		
		elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
			U12_GTAG.append(row)
		
		else:
			indef.append(row)
			
	
	U2 = U2_GTAG
	U12 = U12_ATAC + U12_GTAG
	
	for row in U2:
		dn = row[6]
		dn_U2[dn] += 1
			
	for row in U12:
		dn = row[6]	
		dn_U12[dn] += 1
	
	dn_total_list = dn_total.items()
	dn_U2_list = dn_U2.items()
	dn_U12_list = dn_U12.items()
	
	dn_total_list.sort(key=lambda x: x[1])
	dn_U2_list.sort(key=lambda x: x[1])
	dn_U12_list.sort(key=lambda x: x[1])
	
	print "TOTAL"
	
	for i in reversed(dn_total_list):
		dn = i[0]
		N = i[1]
		print dn, N, percent(N, len(U2_GTAG + U12_ATAC + U12_GTAG + indef))

	print "U2"
	
	for i in reversed(dn_U2_list):
		dn = i[0]
		N = i[1]
		print dn, N, percent(N, len(U2))
		
	print "U12"
	
	for i in reversed(dn_U12_list):
		dn = i[0]
		N = i[1]
		print dn, N, percent(N, len(U12))										
	
	   						

if __name__ == '__main__':
	main(sys.argv[1])	
