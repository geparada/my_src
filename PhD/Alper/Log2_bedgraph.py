import sys
import csv
import numpy as np
import pandas as pd
from collections import defaultdict


def genome_pseudo_counts(chrom_sizes):

	pseudo_counts = defaultdict(int)

	for row in csv.reader(open(chrom_sizes), delimiter = '\t'):

		chrom, size = row

		size = int(size)

		for i in range(size):
			try:

				start = i
				end = i + 1

				index = "-".join(map(str, [chrom, start, end]))
				pseudo_counts[index] += 1

			except TypeError:  # for some reazon i have a problem to process the row after ['chrIII', '13301990', '13301991']
				pass

	return pseudo_counts



def bedgraph_indexer(bedgraph, pseudo_counts):

	counts = pseudo_counts

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



def main(bedgraph_wt, bedgraph_mut, chrom_sizes):

	pseudo_counts_wt = genome_pseudo_counts(chrom_sizes)
	pseudo_counts_mut = pseudo_counts_wt.copy()

	counts_wt, index_wt = bedgraph_indexer(bedgraph_wt, pseudo_counts_wt)
	counts_mut, index_mut = bedgraph_indexer(bedgraph_mut, pseudo_counts_mut)


	bedgraphs = {"WT" : pd.Series(counts_wt, index=index_wt), "MUT" : pd.Series(counts_mut, index=index_mut)}


	bedgraph_df = pd.DataFrame(bedgraphs)
	bedgraph_df = bedgraph_df.dropna() #Remove NaA
	sum_wt = bedgraph_df.sum(axis=0)["WT"]
	sum_mut = bedgraph_df.sum(axis=0)["MUT"]
	
	bedgraph_df['log2'] = np.log2(( (bedgraph_df["MUT"]/sum_mut)  /  (bedgraph_df["WT"]/sum_wt) ))

	#bedgraph_df['log2'] = np.log2(( (bedgraph_df["MUT"])  /  (bedgraph_df["WT"]) ))


	for index, row in bedgraph_df.iterrows():

		chrom, start, end = index.split("-")
		log2 = row['log2']
		out = map(str, [chrom, start, end, log2])
		print  "\t".join(out)


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])


