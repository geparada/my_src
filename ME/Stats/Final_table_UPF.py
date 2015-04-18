import sys
import csv
import numpy as np





def main(Final_ME_TABLE):

	for row in csv.reader(open(Final_ME_TABLE), delimiter = ' '):

		ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, ME_cov_sum, SJ_cov_sum,  mixture_cov, adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov, HepG2_control_cov, HepG2_UPF2_cov, HELA_control_cov, HELA_UPF1_cov = row

		HepG2_ME_control, HepG2_SJ_control, HepG2_ME_UPF, HepG2_SJ_UPF = map(int, HepG2_control_cov.split("/")) + map(int, HepG2_UPF2_cov.split("/"))

		HepG2_control_total = HepG2_ME_control + HepG2_SJ_control
		HepG2_UPF_total =  HepG2_ME_UPF + HepG2_SJ_UPF

		if HepG2_control_total>10 and HepG2_UPF_total>10:

			HepG2_control_phi = float(HepG2_ME_control) / float(HepG2_control_total)
			HepG2_UPF_phi = float(HepG2_ME_UPF) / float(HepG2_UPF_total)

			if HepG2_control_phi!=0 and  HepG2_UPF_phi!=0:

				print HepG2_control_phi, HepG2_UPF_phi, HepG2_UPF_phi/HepG2_control_phi, np.log2(HepG2_UPF_phi/HepG2_control_phi)






if __name__ == '__main__':
	main(sys.argv[1])
