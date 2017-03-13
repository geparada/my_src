import sys
import csv
from collections import defaultdict
import numpy as np



def main(SAM_introns):


	#introns = defaultdict(set)

	total_cov = 0

	for row in csv.reader(open(SAM_introns), delimiter = ' '):


		chrom, start, end, cov = row

		total_cov += 1


	for row in csv.reader(open(SAM_introns), delimiter = ' '):

		chrom, start, end, cov = row
		start= int(start)
		end = int(end)

		norm_cov = (10**6)*float(cov)/float(total_cov)

		if int(cov)>=3:

			print "\t".join(map(str, (chrom, start-8, end+8, round(norm_cov, 2), 0, ".", start-8, end+8, "0,0,0", 2, "8,8", ",".join(map(str, [0, end-start+8])))))



if __name__ == '__main__':
	main(sys.argv[1])	