import sys
import csv
from collections import defaultdict 


def main(gencode_gff, SS_count):


	SS_counts = defaultdict(int)


	for row in csv.reader(open(SS_count), delimiter = ' '):


		SS, count = row
		count = int(count)

		SS_counts[SS] += count



	for row in csv.reader(open(gencode_gff), delimiter = '\t'):

		if row[0][0]!='#':

			chrom, gff_file, feature, start, end, dot1, strand, dot2, IDs = row

			chrom = chrom.strip("chr")

			if feature == "exon":


				print feature, chrom, start, end, SS_count["_".join([chrom, start])], SS_count["_".join([chrom, end])]



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])

#chr1    HAVANA  exon    3205901 3207317 .       -       .       gene_id "ENSMUSG00000051951.5"; transcript_id "ENSMUST00000162897.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "Xkr4"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "Xkr4-003"; exon_number 2; exon_id "ENSMUSE00000866652.1"; level 2; transcript_support_level "1"; havana_gene "OTTMUSG00000026353.2"; havana_transcript "OTTMUST00000086625.1";