import sys
import csv
from collections import defaultdict 


def main(gencode_gff, SS_count):


	SS_counts = defaultdict(int)   #I choose default dict in this case to give 0 as default value.


	for row in csv.reader(open(SS_count), delimiter = ' '):


		SS, count = row
		count = int(count)

		SS_counts[SS] += count


	transcripts = defaultdict(list)
	exons_5 = defaultdict(set)
	exons_3 =  defaultdict(set)

	for row in csv.reader(open(gencode_gff), delimiter = '\t'):

		if row[0][0]!='#':

			chrom, gff_file, feature, start, end, dot1, strand, dot2, IDs = row

			chrom = chrom.strip("chr")

			if feature == "transcript":

				transcript = " ".join(row)

			if feature == "exon":

				estart = "_".join([chrom, str(int(start)-1)])
				eend = "_".join([chrom, end])

				exons_5[estart].add((eend,  SS_counts[eend]))
				exons_3[eend].add((estart, SS_counts[estart]))

				transcript = "\t".join(row)
				exon = "\t".join(row)

				transcripts[transcript].append(exon)


	max_e_5s = {}
	max_e_3s =  {}

	for i in exons_5.items():

		e_5, e_3_counts = i
		e_3 =  max(e_3_counts,key=lambda item:item[1])
		max_e_5s[e_5] = e_3_counts

	for i in exons_3.items():

		e_3, e_5_counts = i
		e_5 =  max(e_5_counts,key=lambda item:item[1])
		max_e_3s[e_3] = e_5_counts


	for t in transcripts.items():
		transcritp,  exons =  t

		filtered_exons = []

		for e in exons:


			chrom, gff_file, feature, start, end, dot1, strand, dot2, IDs = e.split("\t")

			estart = "_".join([chrom, str(int(start)-1)])
			eend = "_".join([chrom, end])


			max_e_3, max_e_3_count = max_e_5s[estart]
			max_e_5, max_e_5_count = max_e_3s[eend]


			if  max_e_5s[estart]==eend and max_e_3s[eend]==estart and  max_e_5_count > 0 and max_e_3_count > 0:

				filtered_exons.append(e)

		if len(filtered_exons)>=0:

			print transcritp

			for e in filtered_exons:

				print e



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])

#chr1    HAVANA  exon    3205901 3207317 .       -       .       gene_id "ENSMUSG00000051951.5"; transcript_id "ENSMUST00000162897.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "Xkr4"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "Xkr4-003"; exon_number 2; exon_id "ENSMUSE00000866652.1"; level 2; transcript_support_level "1"; havana_gene "OTTMUSG00000026353.2"; havana_transcript "OTTMUST00000086625.1";