import sys
import csv


def main(sim_fastq):

	MEs = set([])

	for row in csv.reader(open(sim_fastq), delimiter = '\t'):


		if row[0][1]=="@":

			SJ, ME_seq, ME_start, ME_end, total_coverage, n = row.split("_")

			chr = SJ.split(":")[1:]

			ME_start = int(ME_end)
			ME_end = int(ME_end)

			MEs.add((chr, ME_start, ME_end))


	print MEs










if __name__ == '__main__':
	main (sys.argv[1]) 		