import sys
import csv


def main(DEXSeq_out, gff_chunk):

	chunk_exons = []

	for row in csv.reader(open(gff_chunk), delimiter = '\t'):

		chrom, gff_file, feature, start, end, dot1, strand, dot2, IDs = row

		if feature=="exonic_part":

			transcripts, transcripts_val, exonic_part_number, exonic_part_number_val, gene_id, gene_id_val =  IDs.split(" ")

			chunk_exons.append(":".join([gene_id_val.strip(';"'), exonic_part_number_val.strip(';"')]) )



	exon_count_dict = {}

	for row in csv.reader(open(DEXSeq_out), delimiter = '\t'):


		exon, count = row



		exon_count_dict[exon] = count



	for exon in chunk_exons:

		# try:

		print exon, exon_count_dict[exon]

		# except KeyError:

		# 	pass








if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])