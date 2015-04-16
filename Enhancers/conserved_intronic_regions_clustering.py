import sys
import csv
from collections import defaultdict

def main(TF_cluster):

	clusters = defaultdict(list)

	for row in csv.reader(open(TF_cluster), delimiter = '\t'):

		#print row

		chrom, chromStart, chromEnd, name, score, strand, cluster_ID = row

		chromStart = int(chromStart)
		chromEnd = int(chromEnd)

		clusters[cluster_ID].append((chrom, chromStart, chromEnd, name))


	for k  in clusters.items():

		names = set([])

		chromStarts = []
		chromEnds = []
		chroms = set([])

		cluster_ID, TFs = k

		for i in TFs:

			chrom, chromStart, chromEnd, name = i

			names.add(name.split("|")[0].split("-")[0])
			chromStarts.append(chromStart)
			chromEnds.append(chromEnd)
			chroms.add(chrom)

		chromStart = min(chromStarts)
		chromEnd = max(chromEnds)
		chroms = str(chroms)
		name = cluster_ID + "_" + ",".join(names)


		BED = [chrom, chromStart, chromEnd, name, chromEnd-chromStart, "."]

		print "\t".join(map(str, BED))


 




if __name__ == '__main__':
	main(sys.argv[1])