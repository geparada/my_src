import sys
import csv
import wWigIO
from collections import defaultdict
from ngslib import BigWigFile




def bigwig_fetcher(bigwig, ichr, istart, iend):

	bw=BigWigFile(bigwig)
	
	scores = []

	with BigWigFile(bigwig) as bw:

		for i in bw.fetch(chrom=ichr,start=istart,stop=iend):

			scores.append(i.score)

		return scores

	bw.wWigIO.close()	
	



def main(gencode_bed, phylop_vertebrates):
	
	min_phylop_score = 4.88
	seed_lenght = 50

	for row in csv.reader(open(gencode_bed), delimiter = '\t'):
		
		csv.field_size_limit(1000000000)

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		ichr = row[0]
		gene_name = row[6]

		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = ichr + ":" +  str(istart) + strand + str(iend)

			intron_scores = bigwig_fetcher(phylop_vertebrates, ichr, istart, iend)

			C = 15

			while C <= (iend -istart - 15 - seed_lenght):

				mean_conservation_seed = sum(intron_scores[C:C+seed_lenght])/seed_lenght

				#print mean_conservation_seed, len(intron_scores), C

				if mean_conservation_seed >= min_phylop_score:

					region_id = gene_name + "|" + ichr + ":" + str(istart + C) + "-" + str(istart + seed_lenght +C) 

					print "\t".join(map(str, ([ichr, istart + C, istart + seed_lenght +C, region_id, mean_conservation_seed, strand]) ))	

				C += 10

				if C > (iend -istart - 30 - seed_lenght):
					break


if __name__ == '__main__':	
	main(sys.argv[1], sys.argv[2])


	#python ~/my_src/Enhancers/conserved_intronic_regions_analysis.py ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 ~/db/hg19.100way.phyloP100way.bw