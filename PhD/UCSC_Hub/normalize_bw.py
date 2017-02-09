import sys
import csv

def main(bedGraph):


	print sum(int(row[3]) for row in csv.reader(open(bedGraph), delimiter = '\t'))

	# for row in csv.reader(open(bedGraph), delimiter = '\t'):

	# 	chrom, start, end = row





if __name__ == '__main__':
	main(sys.argv[1])	