import sys
import csv
from collections import defaultdict
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import wWigIO


def main(row_ME):

	ME_reads = defaultdict(set)
	ME_SJs = defaultdict(set)
	SJ_info = {}

	SJ_SJ = defaultdict(set)
	SJ_same_ME = set([])

	for row in csv.reader(open(row_ME), delimiter = ' '):

		read, seq, qual, tag_alingment, t_score, genome_alingment, g_score, same_ME, len_micro_exon_seq_found, micro_exon_seq_found, number_of_micro_exons_matches, max_U2_scores, max_mean_conservations_vertebrates, max_mean_conservations_primates, micro_exons_coords, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = row
		SJ =  tag_alingment.split("|")[0]
		info = " ".join([SJ, len_micro_exon_seq_found, micro_exon_seq_found, number_of_micro_exons_matches, max_U2_scores, max_mean_conservations_vertebrates, max_mean_conservations_primates, micro_exons_coords, U2_scores, mean_conservations_vertebrates, mean_conservations_primates])

		ME_reads[info].add(seq)

		for ME in  micro_exons_coords.split(","):
			ME_SJs[ME].add(SJ + "_" + micro_exon_seq_found)
	

###Cheking if any ME is in two or more SJ

	for i in ME_SJs.items():
		ME, SJs = i

		for SJ_A in SJs:

			for SJ_B in SJs:

				SJ_SJ[SJ_A].add(SJ_B)


	for i in SJ_SJ.items():

		SJ_same_ME.add(" ".join(sorted(list(i[1]))))

####Coverage dict #####

	for i in ME_reads.items():

		SJ, len_micro_exon_seq_found, micro_exon_seq_found, number_of_micro_exons_matches, max_U2_scores, max_mean_conservations_vertebrates, max_mean_conservations_primates, micro_exons_coords, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = i[0].split(" ")

		coverage = len(i[1])

		info = " ".join([str(coverage), len_micro_exon_seq_found, micro_exon_seq_found, number_of_micro_exons_matches, max_U2_scores, max_mean_conservations_vertebrates, max_mean_conservations_primates, micro_exons_coords, U2_scores, mean_conservations_vertebrates, mean_conservations_primates])
		SJ_MEseq = "_".join([SJ, micro_exon_seq_found])

		SJ_info[SJ_MEseq] = info


####

	for i in SJ_same_ME:

		#if len(i.split(" ")) > 1:

		SJs = []
		SJ_Coverages = []
		SJ_number_of_micro_exons_matches = []
		SJ_max_U2_scores = []
		SJ_max_mean_conservations_vertebrates = []
		SJ_max_mean_conservations_primates = []
		ME = []

		
		P_MEs  = []

		info = set([])

		for SJ_MEseq in i.split(" "):

			SJ = SJ_MEseq.split("_")[0]

			coverage, len_micro_exon_seq_found, micro_exon_seq_found, number_of_micro_exons_matches, max_U2_scores, max_mean_conservations_vertebrates, max_mean_conservations_primates, micro_exons_coords, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = SJ_info[SJ_MEseq].split(" ")

			SJs.append(SJ)
			SJ_Coverages.append(int(coverage))
			info.add((len_micro_exon_seq_found, micro_exon_seq_found))
			SJ_number_of_micro_exons_matches.append(int(number_of_micro_exons_matches))
			SJ_max_U2_scores.append(float(max_U2_scores))
			SJ_max_mean_conservations_vertebrates.append(float(max_mean_conservations_vertebrates))
			SJ_max_mean_conservations_primates.append(float(max_mean_conservations_primates))

			SJ_chr, SJ_istart, SJ_iend = re.findall(r"[\w']+", SJ)
			SJ_istart = int(SJ_istart)
			SJ_iend = int(SJ_iend)

			len_micro_exon_seq_found = int(len_micro_exon_seq_found)

			SJ_len = SJ_iend - SJ_istart
			Kmer = SJ_len - (len_micro_exon_seq_found+1)
			P_ME = 1 - ( 1 - (float(1)/float(4**len_micro_exon_seq_found+4)))**Kmer	

			P_MEs.append(P_ME)

			set_ME = set([])

			for a, b, c, d in zip(micro_exons_coords.split(","), U2_scores.split(","), mean_conservations_vertebrates.split(","), mean_conservations_primates.split(",")):

				set_ME.add("|".join([a,b,c, d]))

			ME.append(set_ME)

		sum_total_coverage = sum(SJ_Coverages)
		total_SJs = ",".join(SJs)
		total_coverages = ",".join(map(str, SJ_Coverages))
		
		total_max_U2_scores = min(SJ_max_U2_scores)
		total_max_mean_conservations_vertebrates = min(SJ_max_mean_conservations_vertebrates)
		total_max_mean_conservations_primates = min(SJ_max_mean_conservations_primates)


		total_ME = ",".join(set.intersection(*ME))
		total_number_of_micro_exons_matches = len(total_ME.split(","))

		len_micro_exon_seq_found, micro_exon_seq_found = list(info)[0]

		#if 6 >= len(micro_exon_seq_found) >= 3:

		#### Probabilidad ###


		out =  map(str, [sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min(P_MEs), total_ME])
		
		if total_ME=="":

			print "\t".join(out)





if __name__ == '__main__':	
	main(sys.argv[1])


#python ~/my_src/ME/Pipeline/ME_centric_table.py _clip1.trim.sam.row_ME.filter1