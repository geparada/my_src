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
		strand = row[5]
		EST_coverage = row[7]
		cDNA_coverage = row[8]
		EST = row[9]
		cDNA = row[10]
		intron = chr + ":" + istart + strand + iend
		TOTAL_coverage =  int(EST_coverage) + int(cDNA_coverage)
		
		if TOTAL_coverage >= 3:
			print chr, istart, iend, "".join(( intron, "|", dn,  "|", cDNA_coverage, "|", EST_coverage, "|", cDNA, "|", EST)), 0, strand    
	
	
	
if __name__ == '__main__':
	main(sys.argv[1])
