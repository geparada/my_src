import sys
import csv
import wWigIO
from collections import defaultdict
from ngslib import BigWigFile


phylop_bigwig = defaultdict(list)

def bigwig_loader(bigwig, chromInfo):

	print >> sys.stderr, "Loading the BigWig file in RAM memory ...",		

	bw=BigWigFile(bigwig)
	
	for row in csv.reader(open(chromInfo), delimiter = '\t'):

		chr = row[0]
		chr_size = int(row[1])

		for i in bw.fetch(chrom=chr,start=0,stop=chr_size):

			phylop_bigwig[chr].append(i.score)

		bw.wWigIO.close()

	print >> sys.stderr, "OK",	


def bigwig_mean(bigwig, chr, start, end):


	score_sum = 0
	mean_score = 0

	if (end-start) != 0:

		mean_score = sum(bigwig[chr][start:end])/(end-start)

	return mean_score



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
			intron = ichr + ":" +  str(istart) + row[5] + str(iend)

			C = 0			

			while C <= (iend -istart - 30 - seed_lenght):

				mean_conservation_seed = bigwig_mean(phylop_bigwig, ichr, C+istart, C+istart+seed_lenght)

				if mean_conservation_seed >= min_phylop_score:

					region_id = gene_name + "|" + ichr + ":" + str(istart + C) + "-" + str(istart + seed_lenght +C) 

					print "\t".join(map(str, ([ichr, istart + C, istart + seed_lenght +C, region_id, mean_conservation_seed, strand]) ))	

				C += 10

				if C > (iend -istart - 30 - seed_lenght):
					break


if __name__ == '__main__':
	bigwig_loader(sys.argv[2],sys.argv[3])	
	main(sys.argv[1])


	#python ~/my_src/Enhancers/conserved_intronic_regions_analysis.py ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 ~/db/hg19.100way.phyloP100way.bw