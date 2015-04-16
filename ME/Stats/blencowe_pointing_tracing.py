import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


     

def main(filter3, filter2, filter1, out_filter2, row_list, blencowe, ponting):

	ME_row = set([])
	ME_filter1 = set([])
	ME_filter2 = set([])
	ME_filter3 = set([])

	SNP = set([])
	homopolymer = set([])
	Burge = set([])
	low_coverage = set([])

	blencowe_row = set([])
	blencowe_filter1 = set([])
	blencowe_filter2 = set([])
	blencowe_filter3 = set([])
	blencowe_SNP = set([])
	blencowe_homopolymer = set([])
	blencowe_Burge = set([])
	blencowe_low_coverage = set([])

	ponting_row = set([])
	ponting_filter1 = set([])
	ponting_filter2 = set([])
	ponting_filter3 = set([])
	ponting_SNP = set([])
	ponting_homopolymer = set([])
	ponting_Burge = set([])
	ponting_low_coverage = set([])

	for row in csv.reader(open(row_list), delimiter = '\t'):

		read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found, micro_exons_coords = row

		for ME in micro_exons_coords.split(","):

			ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")
			ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end
			ME_row.add(ME_ID)


	for row in csv.reader(open(filter1), delimiter = ' '):

		read, seq, qual, tag_alingment, t_score, genome_alingment, g_score, same_ME, len_DR_corrected_micro_exon_seq_found, DR_corrected_micro_exon_seq_found, len_micro_exons, max_U2_scores, max_TOTAL_mean_conservation_vertebrates, max_TOTAL_mean_conservation_primates, micro_exons_coords, U2_scores, TOTAL_mean_conservation_vertebrates, TOTAL_mean_conservation_primates = row

		for ME in micro_exons_coords.split(","):

			ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")
			ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end
			ME_filter1.add(ME_ID)

	for row in csv.reader(open(out_filter2), delimiter = '\t'):

		#read, seq, qual, tag_alingment, t_score, genome_alingment, g_score, same_ME, len_DR_corrected_micro_exon_seq_found, DR_corrected_micro_exon_seq_found, len_micro_exons, max_U2_scores, max_TOTAL_mean_conservation_vertebrates, max_TOTAL_mean_conservation_primates, micro_exons_coords, U2_scores, TOTAL_mean_conservation_vertebrates, TOTAL_mean_conservation_primates, out_reason = row

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, P_ME, total_ME, out_reason = row

		#for ME in micro_exons_coords.split(","):

		for ME_scores in total_ME.split(","):

			ME = ME_scores.split("|")[0]

			if ME!="":

				ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")
				ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end

				if out_reason=="SNP":

					SNP.add(ME_ID)

				if out_reason=="homopolymer":
					homopolymer.add(ME_ID)

				if out_reason=="Burge":
					Burge.add(ME_ID)

				if out_reason=="coverage":
					low_coverage.add(ME_ID)


	for row in csv.reader(open(filter2), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, P_ME, total_ME = row

		for ME_scores in total_ME.split(","):

			ME = ME_scores.split("|")[0]
			#print ME
			if ME!="":

				ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")
				ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end
				ME_filter2.add(ME_ID)


	for row in csv.reader(open(filter3), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_ME, total_ME, true_ME, score, is_annotated = row

		ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = true_ME.split("|")
		ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")
		ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end
		ME_filter3.add(ME_ID)

	for row in csv.reader(open(blencowe), delimiter = '\t'):

		if int(row[3])<=25:

			ME_chrom = row[2].split(":")[0]
			ME_start = str(int(row[2].split(":")[1].split("-")[0]) -1)
			ME_end = row[2].split(":")[1].split("-")[1]
			ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end

			if ME_ID in ME_filter3:
				blencowe_filter3.add(ME_ID)

			elif ME_ID in ME_filter2:
				blencowe_filter2.add(ME_ID)

			elif ME_ID in ME_filter1:
				blencowe_filter1.add(ME_ID)

				if ME_ID in SNP:
					blencowe_SNP.add(ME_ID)
				if ME_ID in homopolymer:
					blencowe_homopolymer.add(ME_ID)
				if ME_ID in Burge:
					blencowe_Burge.add(ME_ID)
				if ME_ID in low_coverage:
					blencowe_low_coverage.add(ME_ID)

			elif ME_ID in ME_row:
				blencowe_row.add(ME_ID)


	for row in csv.reader(open(ponting), delimiter = '\t'):

		if int(row[4])<=25:

			ME_chrom = row[0]
			ME_start = row[1]
			ME_end = row[2]
			ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end

			if ME_ID in ME_filter3:
				ponting_filter3.add(ME_ID)

			elif ME_ID in ME_filter2:
				ponting_filter2.add(ME_ID)

			elif ME_ID in ME_filter1:
				ponting_filter1.add(ME_ID)

				if ME_ID in SNP:
					ponting_SNP.add(ME_ID)
				if ME_ID in homopolymer:
					ponting_homopolymer.add(ME_ID)
				if ME_ID in Burge:
					ponting_Burge.add(ME_ID)
				if ME_ID in low_coverage:
					ponting_low_coverage.add(ME_ID)

			elif ME_ID in ME_row:
				ponting_row.add(ME_ID)

	print "List", "Blencowe_et_al", "Ponting_et_al"

	print "Final_List", len(blencowe_filter3), len(ponting_filter3)

	print "Pre_Filter_3", len(blencowe_filter2), len(ponting_filter2)

	print "Pre_filter_2", len(blencowe_filter1), len(ponting_filter1)

	print "Pre_filter_1", len(blencowe_row), len(ponting_row)

	print "\n"

	print "Filter_2_out_reason", "Blencowe_et_al", "Ponting_et_al"

	print "SNP", len(blencowe_SNP), len(ponting_SNP)

	print "Homopolymer", len(blencowe_homopolymer), len(ponting_homopolymer)

	print "Burge", len(blencowe_Burge), len(ponting_Burge)

	print "Low_coverage", len(blencowe_low_coverage), len(ponting_low_coverage)



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])

#python ~/my_src/ME/Stats/blencowe_pointing_tracing.py TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 TOTAL.sam.row_ME.filter1.ME_centric.filter2 TOTAL.sam.row_ME.filter1 TOTAL.sam.row_ME.filter1.ME_centric.filter2.out TOTAL.sam.row_ME ../../Literature/Blencowe_et_al.txt.csv ../../Literature/Ponting_et_al.txt