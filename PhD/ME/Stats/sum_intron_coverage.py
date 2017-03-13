import sys
import csv
from collections import defaultdict
import numpy as np




intron_coverage = defaultdict(int)

for i in range(1, len(sys.argv)):

	for row in csv.reader(open(sys.argv[i]), delimiter = ' '):

		chrom, start, end, cov = row
		intron_coverage[" ".join(row[:-1])] += int(cov)


for i in intron_coverage.items():

	print i[0], i[1]