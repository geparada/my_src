import sys
import csv
from collections import defaultdict

def main(ME_centric):

	ME_lens = defaultdict(int)
	

	ME_len_cov = defaultdict(list)

	for row in csv.reader(open(ME_centric), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, total_ME = row

		len_micro_exon_seq_found = int(len_micro_exon_seq_found)
		total_max_U2_scores = float(total_max_U2_scores)
		sum_total_coverage = int(sum_total_coverage)

		# ME_lens[len_micro_exon_seq_found] += 1
		# ME_cov[sum_total_coverage] += 1

		ME_len_cov[len_micro_exon_seq_found].append(sum_total_coverage)


	for i in ME_len_cov.items():

		ME_len, ME_cov = i

		S_0_5 = 0
		S_5_10 = 0
		S_10_20 = 0
		S_20_40 = 0
		S_40_ = 0

		for c in ME_cov:

			if 0 < c <=5:
				S_0_5 += 1

			elif 5 < c <=10:
				S_5_10 += 1

			elif 10 < c <=20:
				S_10_20 += 1

			elif 20 < c <=40:
				S_20_40 += 1

			elif 40 < c:
				S_40_ += 1

		print ME_len, S_0_5, S_5_10, S_10_20, S_20_40, S_40_





if __name__ == '__main__':
	main(sys.argv[1])	
