import sys
import csv

def main (cDNA_EST, FINAL_TABLE):
	""" Extrae cDNAs/ESTs que apoyan la conservacion de los intrones no canonicos """

	csv.field_size_limit(1000000000)

	non_canonical_introns = set([])


	for row in csv.reader(open(FINAL_TABLE), delimiter = '\t');  
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
		tissues_coverage = int(row[16])
		n_tissues = row[17]
		tissues = row[18]

		if dn !="GTAG" and dn !="GCAG" and dn !="ATAC":
			non_canonical_introns.add(intron) 		

		
	for row in csv.reader(open(file), delimiter = ' '):
		
		intron = row[0]
		coverage = int(row[1])
		chr = row[2]
		strand = row[3]
		istart = row[4]
		iend = row[5]
		ilength = int(row[6])
		dn = row[7]
		dn_type = row[8]
		dn_type_score = row[9]
		reads = row[10].split(",")
		

	
	




if __name__ == '__main__':
	main(sys.argv[1])
