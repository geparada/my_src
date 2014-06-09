import sys
import csv

def main (old_FINAL_TABLE, new_FINAL_TABLE):	
	""" Arregla el coverage de la conservacion """

	intron_mm9 = {}

	for row in csv.reader(open(new_FINAL_TABLE), delimiter = '\t'):
		intron = row[0]
		mm9_cDNA_coverage = row[13]
		mm9_EST_coverage = row[14]
		
		intron_mm9[intron] = [mm9_cDNA_coverage, mm9_EST_coverage]

	
	for row in csv.reader(open(old_FINAL_TABLE), delimiter = '\t'):
		intron = row[0]
		new_mm9 = intron_mm9[intron]
		
		print "\t".join(row[:13] + new_mm9 + row[15:])
#		print "\t".join(map(str, (row[:12] + new_mm9 + row[14:])))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
