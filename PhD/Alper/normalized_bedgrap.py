import sys
import csv
import numpy as np
import pandas as pd
from collections import defaultdict




def bedgraph_indexer(bedgraph):

	counts = defaultdict(int)

	count_list = []
	index_list = []




	for row in csv.reader(open(bedgraph), delimiter = '\t'):

		chrom, start, end, count = row

		start = int(start)
		end = int(end)
		count = float(count)
		index = "-".join(map(str, [chrom, start, end]))

		if end - start == 1:  # If one seqment is more than 1 base lenght, this step descompose it into individual bases

			counts[index] += count

		else:

			for i in range(end - start + 1):

				index = "-".join(map(str, [chrom, start + i , start + i + 1]))
				counts[index] += count

	for row in counts.items():

		index, count = row
		count_list.append(count)
		index_list.append(index)

	return count_list, index_list



def main(bedgraph, chrom_sizes):



	counts, index = bedgraph_indexer(bedgraph)

	bedgraphs = {"bg" : pd.Series(counts, index=index)}
	bedgraph_df = pd.DataFrame(bedgraphs)
	sum_bg= bedgraph_df.sum(axis=0)["bg"]
	
	bedgraph_df['count_norm'] = bedgraph_df["bg"]/sum_bg


	black_list = set([])

	for row in csv.reader(open(bedgraph), delimiter = '\t'):

		chrom, size = row
		size = isize
		black_list.add(chrom + "_" + size)		


	for index, row in bedgraph_df.iterrows():

		chrom, start, end = index.split("-")
		count_norm = row['count_norm']
		out = map(str, [chrom, start, end, count_norm])

		if ( (chrom + "_" +str(end)) in black_list) == False:
			print  "\t".join(out)


if __name__ == '__main__':
	main(sys.argv[1])
