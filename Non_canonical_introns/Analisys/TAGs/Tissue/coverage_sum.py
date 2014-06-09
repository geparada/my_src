import sys
import csv
csv.field_size_limit(1000000000)
from collections import defaultdict

intron_coverage = defaultdict(int)

def coverage_sum(SJ_coverage):
	""" Extrae el coverage de los distintos archivos """
	
	reader1 = csv.reader(open(SJ_coverage), delimiter=" ")
	
	for row in reader1:
		intron = row[0]
		coverage = int(row[1])
		intron_coverage[intron] += coverage	
	

def main (s1x75, p2x50_1, p2x50_2):
	""" Suma los coverage dentro de cada tejido """
	
	coverage_sum(s1x75)
	coverage_sum(p2x50_1)	
	coverage_sum(p2x50_2)
	
	for i in intron_coverage.items():
		print i[0], i[1]  	
	

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])
