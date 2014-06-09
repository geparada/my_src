import sys
import csv
from collections import defaultdict
import scipy.stats
import numpy as np
from numpy import array


def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = ' ')

	GTAG_U2 = defaultdict(int)
	ATAC_U12 = defaultdict(int)
	GTAG_U12 = defaultdict(int)
	

	GTAG_U2_list = []
	ATAC_U12_list = []
	GTAG_U12_list = []

	for row in reader1:
		intron = row[0]
		coverage = row[1]
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		ilength = row[6]
		dn = row[7]
		dn_type = row[8]
		dn_type_score = float(row[9])
		
		
		if dn_type == "U2_GTAG":
			GTAG_U2[int(dn_type_score)] += 1
			GTAG_U2_list.append(dn_type_score)

		if dn_type == "U12_ATAC":		
			ATAC_U12[int(dn_type_score)] += 1
			ATAC_U12_list.append(dn_type_score)					
		
		if dn_type == "U12_GTAG":
			GTAG_U12[int(dn_type_score)] += 1
			GTAG_U12_list.append(dn_type_score)
	
	
	print "dn_type", "(media + 2std)", "top5%"


	GTAG_U2_list.sort()
	ATAC_U12_list.sort()
	GTAG_U12_list.sort()	
	GTAG_U2_array = np.array(GTAG_U2_list)
	ATAC_U12_array = np.array(ATAC_U12_list)	
	GTAG_U12_array = np.array(GTAG_U12_list)
	
	print "U2_GTAG", GTAG_U2_array.mean() + 2 * GTAG_U2_array.std(), GTAG_U2_list[int(round(-len(GTAG_U2_list)/20))]
	print "U12_ATAC", ATAC_U12_array.mean() + 2 * ATAC_U12_array.std(), ATAC_U12_list[int(round(-len(ATAC_U12_list)/20))]
	print "U12_GTAG", GTAG_U12_array.mean() + 2 * GTAG_U12_array.std(), GTAG_U12_list[int(round(-len(GTAG_U12_list)/20))]	
	
#	print np.array(GTAG_U2_list).std()
	
	print "dn_type_score", "U2_GTAG", "U12_ATAC", "U12_GTAG"
	
	for i in range(100):
		print i, GTAG_U2[i], ATAC_U12[i], GTAG_U12[i]

#Almacenar los valores en los defaults dicts utilizando los int(dn_type_score) como key

if __name__ == '__main__':
	main(sys.argv[1])
