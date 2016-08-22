import sys
import csv



def main(narrowPeak):


	for i in  in csv.reader(open(narrowPeak), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak = row

		print chrom, chromStart, chromEnd, name, score,



if __name__ == '__main__':
	main(sys.argv[1])

