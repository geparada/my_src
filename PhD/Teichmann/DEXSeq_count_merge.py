import sys
import csv
from collections import defaultdict
import numpy as np



cell_exon_count = defaultdict(list)
cell_gene_count = defaultdict(list)
exon_number_exons_per_gene = defaultdict(set)


if __name__ == '__main__':


	gene_lenghts  = defaultdict(int)

	for row in csv.reader(open(sys.argv[1]), delimiter = '\t'):

		chrom, gff_file, feature, start, end, dot1, strand, dot2, IDs = row

		if feature == "exonic_part":

			start = int(start)
			end = int(end)

			lenght = end - start

			gene = IDs.split(" ")[-1].strip('"')

			gene_lenghts[gene] += lenght



	cells = [i.split(".")[0] for i in sys.argv[2:] ] 


	library_sizes = []

	for i in range(2, len(sys.argv)):

		gene_count = defaultdict(int)

		library_size = 0


		for row in csv.reader(open(sys.argv[i]), delimiter = '\t'):

			exon, count = row

			count = int(count)

			if exon[0]!="_":

				cell_exon_count[exon].append(count)

				gene, exon_number = exon.split(":")

				gene_count[gene]+= count

				exon_number_exons_per_gene[gene].add(exon_number)

				library_size += count


		library_sizes.append(library_size)

		for i in gene_count.items():

			gene, count =  i
			cell_gene_count[gene].append(count)


	for i in cell_exon_count.items():

		

		exon, exon_counts = i

		gene, exon_number = exon.split(":")


		gene_counts =  cell_gene_count[gene]

		ratios = []
		phis = []
		rpkms = []
		

		for e, g, lib_size, cell in zip(exon_counts, gene_counts, library_sizes, cells):

			if e!=0:


				
				ratio = str(e) + "/" + str(g)
				phi= float(e)/float(g)
				gene_lenght = gene_lenghts[gene]


				rpkm = (10**9 * float(g))/ (float(gene_lenght) * float(lib_size))




				#if rpkm >= 10 and e >= 10:
				#if rpkm >= 10:

				ratios.append(ratio)
				phis.append(phi)
				rpkms.append(rpkm)

		CV_phis =  (100* np.std(phis)) / np.mean(phis)

		if len(phis)>5:

			# print cell, exon, gene, exon_number, ratios, phis, rpkms, CV_phis, np.mean(rpkms)

			exon_number = "E" + exon_number

			print exon, gene, exon_number, CV_phis, np.mean(rpkms), ",".join(ratios), ",".join(map(str, phis)), ",".join(map(str, rpkms))


# python ~/my_src/PhD/Teichmann/DEXSeq_count_merge.py gencode.vM9.chr_patch_hapl_scaff.annotation.DEXseq.gff.sed i2m_* | sort > i2m__phi_CV.txt



