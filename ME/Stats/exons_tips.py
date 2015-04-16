import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


     

def main(bed12, ME_centric_filter3, blencowe, ponting):

	estarts = set([])
	eends = set([])
	blencowe_ME = set([])
	ponting_ME = set([])

	for row in csv.reader(open(blencowe), delimiter = '\t'):

		blencowe_ME.add(row[2])


	for row in csv.reader(open(ponting), delimiter = '\t'):

		ME_chrom = row[0]
		ME_start = row[1]
		ME_end = row[2]
		ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end
		ponting_ME.add(ME_ID)


	for row in csv.reader(open(bed12), delimiter = '\t'):
		
		csv.field_size_limit(1000000000)

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]


		for q1, b in zip(qstarts, blocksizes):
			estart = start + q1
			eend = start + q1 + b
			elenght = eend - estart
			exon = chr + ":" +  str(estart) + strand + str(eend)

			if eend - estart > 25:

				if strand == "+":

					estarts.add(chr + "_" +  str(estart))
					eends.add(chr + "_" +  str(eend))

				elif strand == "-":

					eends.add(chr + "_" +  str(estart))
					estarts.add(chr + "_" +  str(eend))


	for row in csv.reader(open(ME_centric_filter3), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_ME, total_ME, true_ME, score, is_annotated = row

		ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = true_ME.split("|")

		ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")

		len_micro_exon_seq_found = int(len_micro_exon_seq_found)

		alternative_5 = False
		alternative_3 = False

		blencowe_annotated = False
		ponting_annotated = False


		if (ME_chrom + ":" + str(int(ME_start)+1) + "-" + ME_end) in blencowe_ME:
			blencowe_annotated = True

		if (ME_chrom + ":" + ME_start + "-" + ME_end) in ponting_ME:
			ponting_annotated = True


		if ME_strand == "+":

			if (ME_chrom + "_" + ME_start) in estarts:

				alternative_5 = True

			if (ME_chrom + "_" + ME_end) in eends:

				alternative_3 = True

		elif ME_strand == "-":

			if (ME_chrom + "_" + ME_end) in estarts:

				alternative_5 = True

			if (ME_chrom + "_" + ME_start) in eends:

				alternative_3 = True


		if 25 >= len_micro_exon_seq_found >=3:

			print row


			# if blencowe_annotated==True:

			# 	print

		 	# if is_annotated=="False":

		 	# 	# print row

		 	# 	if blencowe_annotated==False and ponting_annotated==False:

		 	# 		if alternative_3:

		 	# 			print row

			 		#print "\t".join(row) + "\t" + "Novel"

			#  	else:

			#  		print "\t".join(row) + "\t" + "Annotated"

			# else:

			# 	print "\t".join(row) + "\t" + "Annotated"

		

		# if is_annotated=="True":

		# 	print row

			# if blencowe_annotated:

			# 	print row



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

#python ~/my_src/ME/Stats/exons_tips.py ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 ../../Literature/Blencowe_et_al.txt.csv ../../Literature/Ponting_et_al.txt