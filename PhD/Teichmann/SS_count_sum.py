import sys
import csv
from collections import defaultdict


SS_sum = defaultdict(int)

for i in range(1, len(sys.argv)):

	for row in csv.reader(open(sys.argv[i]), delimiter = ' '):

		SS, count = row
		count = int(count)

		SS_sum[SS] += count


for i in SS_sum.items():

	print " ".join(i)





