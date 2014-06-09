import sys
import csv


def main(TF):
	
	csv.field_size_limit(1000000000)
	reader = csv.reader(open(TF), delimiter = '\t')
	
	for row in reader:
		dn = row[1]
		chr = row[2]
		istart = row[3]
		iend = row[4]
		starnd = row[5]
		EST_coverage = row[7]
		cDNA_coverage = row[8]
		EST = row[9]
		cDNA = row[10]
		
		intron = chr + ':' + istart + starnd + iend
		
		if dn != "GTAG" and dn != "GCAG" and dn != "ATAC":
			print intron, cDNA_coverage, chr, starnd, istart, iend, EST_coverage, dn, cDNA + ',' + EST
		
		else:
			print intron, cDNA_coverage, chr, starnd, istart, iend, EST_coverage, dn, 0
		
if __name__ == '__main__':
	main(sys.argv[1])
