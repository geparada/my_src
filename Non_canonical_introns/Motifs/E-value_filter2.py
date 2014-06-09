import csv
import sys
import re
from collections import defaultdict
import numpy

def main (motifs, sample):
	""" Filtra por E-values los motivos de meme """
	
	reader1 = csv.reader(open(motifs), delimiter = '\t')
	reader2 = csv.reader(open(sample), delimiter = '\t')	

	E_samples = defaultdict(list)

	for row in reader2:
		E = float(row[1].split(" = ")[-1])
		sample = re.findall(r"[\w']+", row[0])[2]
		E_samples[sample].append(E)
	
	min_Es = []
	
	for i in E_samples.items():
		min_Es.append(min(i[1]))
	
	E_threshold = min(min_Es)
	
	for row in reader1:
		E = float(row[1].split(" = ")[-1])
		if E < E_threshold:
			print "\t".join(row)
		


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])	
