import sys
import csv
from scipy import stats

def main(ME_centric, exon_scores):

	gencode_U2_scores = []
	gencode_mean_conservation_vertebrates = []
	gencode_mean_conservation_primates = []


	for row in csv.reader(open(exon_scores), delimiter = ' '):

		chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates = row

		gencode_U2_scores.append(float(U2_score))
		gencode_mean_conservation_vertebrates.append(float(mean_conservation_vertebrates))
		gencode_mean_conservation_primates.append(float(mean_conservation_primates))


	for row in csv.reader(open(ME_centric), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, P_ME, total_ME = row


		total_max_U2_scores = float(total_max_U2_scores)
		total_max_mean_conservations_vertebrates = float(total_max_mean_conservations_vertebrates)
		total_max_mean_conservations_primates = float(total_max_mean_conservations_primates)
		P_ME = float(P_ME)

		percentil_U2_score = stats.percentileofscore(gencode_U2_scores, total_max_U2_scores)
		percentil_mean_conservation_vertebrates = stats.percentileofscore(gencode_mean_conservation_vertebrates, total_max_mean_conservations_vertebrates)
		percentil_mean_conservation_primates = stats.percentileofscore(gencode_mean_conservation_primates, total_max_mean_conservations_primates)


		overall_scores = []

		for m in total_ME.split(","):
			ME, U2_score, mean_conservation_primates, mean_conservation_vertebrates = m.split("|")

			U2_score = float(U2_score)
			mean_conservation_primates = float(mean_conservation_primates)
			mean_conservation_vertebrates = float(mean_conservation_vertebrates)
			
			ME_percentil_U2_score = stats.percentileofscore(gencode_U2_scores, U2_score)
			ME_percentil_mean_conservation_vertebrates = stats.percentileofscore(gencode_mean_conservation_vertebrates, mean_conservation_primates)
			ME_percentil_mean_conservation_primates = stats.percentileofscore(gencode_mean_conservation_primates, mean_conservation_vertebrates)



			overall_score = P_ME * (1- ME_percentil_U2_score/100) * (1 - ME_percentil_mean_conservation_vertebrates/100)

			if ME_percentil_mean_conservation_primates > ME_percentil_mean_conservation_vertebrates:
				overall_score = P_ME * (1- ME_percentil_U2_score/100) * (1 - ME_percentil_mean_conservation_primates/100)


			overall_scores.append(overall_score)

		print "\t".join(map( str, [sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, percentil_U2_score, total_max_mean_conservations_vertebrates, percentil_mean_conservation_vertebrates, total_max_mean_conservations_primates, percentil_mean_conservation_primates, P_ME, min(overall_scores), total_ME ]))

		#"\t".join(map( str, [min(overall_scores), sum_total_coverage, total_SJs, micro_exon_seq_found, total_max_U2_scores, percentil_U2_score, total_max_mean_conservations_vertebrates, percentil_mean_conservation_vertebrates, total_max_mean_conservations_primates, percentil_mean_conservation_primates, P_ME]))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])	