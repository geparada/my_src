import sys
import csv
from collections import defaultdict



def main(gff_dexseq):
	
	chunk_size = 4645
	chunks_row = defaultdict(list)
	g = 0
	chunk = 1


	for row in csv.reader(open(gff_dexseq), delimiter = '\t'):

		chrom, gff_file, feature, start, end, dot1, strand, dot2, IDs = row

		if feature=="aggregate_gene":
			g += 1

		if g>chunk_size:
			chunk += 1
			g = 0

		chunks_row[chunk].append(row)

	for i in chunks_row.items():

		chunk, rows = i
		out = "gff_chunks/" + sys.argv[1] + "." + str(chunk) 
		f = open(out ,'w')

		for row in rows:

			f.write("\t".join(row)+ "\n")





if __name__ == '__main__':
	main(sys.argv[1])