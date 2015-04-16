import sys
import csv
import numpy as np

def log2 (n):

	if n!=0:
		return np.log(n)

	else:
		return n

def ratio(n,d):

	n = float(n)
	d = float(d)

	if d!=0:
		return n / d

	else:
		return n


def main(sh, control):
	
	for row1, row2 in zip( csv.reader(open(sh), delimiter = '\t'), csv.reader(open(control), delimiter = '\t') ):

		ME_ID, ME_strand, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, is_annotated, total_SJs, ME_SJ_coverages_sh, sum_ME_coverage_sh, SJ_coverages_sh, sum_SJ_coverage_sh  = row1

		ME_ID, ME_strand, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, is_annotated, total_SJs, ME_SJ_coverages_control, sum_ME_coverage_control, SJ_coverages_control, sum_SJ_coverage_control  = row2

		sh_phi = ratio(sum_ME_coverage_sh, sum_SJ_coverage_sh)
		control_phi = ratio(sum_ME_coverage_control, sum_SJ_coverage_control)

		change = ratio(sh_phi, control_phi)

		if sh_phi!=0 and control_phi!=0:	

			#print ME_ID, len_micro_exon_seq_found, sum_ME_coverage_sh, sum_SJ_coverage_sh, sum_ME_coverage_control, sum_SJ_coverage_control 
			print ME_ID, len_micro_exon_seq_found, sh_phi, log2(sh_phi), control_phi, log2(control_phi), change, log2(change)

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2] )