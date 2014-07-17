import sys
import csv


def main(FANTOM_Robust_TSS_BED, chrom_sizes):

	chrom_lens = {}

	for row in csv.reader(open(chrom_sizes), delimiter = '\t'):

		chrom, chrom_len = row
		chrom_len = int(chrom_len)

		chrom_lens[chrom] = chrom_len

	for row in csv.reader(open(FANTOM_Robust_TSS_BED), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand = row

		chromStart = int(chromStart)
		chromEnd = int(chromEnd)

		chrom_len = chrom_lens[chrom]

		up = 2000
		down = 500

		new_chromStart = 0
		new_chromEnd = 0

		if strand == "+":

			new_chromStart = chromStart - up
			new_chromEnd = chromEnd + down

		if strand == "-":

			new_chromStart = chromStart - down
			new_chromEnd = chromEnd + up

		if new_chromStart < 0:
			new_chromStart = 0

		if chrom_len < new_chromEnd:
			new_chromEnd = chrom_len 

		BED = [chrom, new_chromStart, new_chromEnd, name, score, strand]

		print "\t".join(map(str, BED))



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])