import sys
import csv

def main(bedGraph):


	total_count = sum(int(row[3]) for row in csv.reader(open(bedGraph), delimiter = '\t'))

	for row in csv.reader(open(bedGraph), delimiter = '\t'):

	 	chrom, start, end, count = row

	 	count = int(count)

	 	normalised_count = str((count/total_count)*10**6)

	 	print "\t".join(chrom, start, end, normalised_count)



if __name__ == '__main__':
	main(sys.argv[1])	