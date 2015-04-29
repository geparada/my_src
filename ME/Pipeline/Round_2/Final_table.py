import csv
import sys

ME_alt = {}

def ME_SJ_coverage_reader(ME_SJ, cov):

	for row in csv.reader(open(ME_SJ), delimiter = '\t'):

		#ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, is_annotated, total_SJs, ME_SJ_coverages, sum_ME_coverage, SJ_coverages, sum_SJ_coverage = row

		ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, is_annotated, total_SJs, ME_SJ_coverages, sum_ME_coverage, SJ_coverages, sum_SJ_coverage, is_alternative_5, is_alternative_3, alternatives_5, cov_alternatives_5, total_cov_alternatives_5, alternatives_3, cov_alternatives_3,  total_cov_alternatives_3 = row



		#psi = float(sum_ME_coverage) / (float(sum_ME_coverage) + float(sum_SJ_coverage))

		ME_alt[ME] = [is_alternative_5, is_alternative_3]

		if is_alternative_5=="True" or is_alternative_3=="True":

			cov[ME] = sum_ME_coverage + "/" + sum_SJ_coverage + "|" + total_cov_alternatives_5 + "|" + total_cov_alternatives_3 #+ "_" + str(psi)

		else:
			cov[ME] = sum_ME_coverage + "/" + sum_SJ_coverage

		#(ME_SJ_coverages, sum_ME_coverage , SJ_coverages, sum_SJ_coverage)

def main( ME_centric_filter3, blencowe, ponting, mixture, adipose, adrenal, brain, breast, colon, heart, kidney, liver, lung, lymph_node, ovary, prostate, skeletal_muscle, testes, thyroid,  white_blood_cells, HepG2_control, HepG2_UPF2, HELA_control, HELA_UPF1):

	blencowe_ME = set([])
	ponting_ME = set([])

	for row in csv.reader(open(blencowe), delimiter = '\t'):
		if int(row[3])<=25:

			ME_chrom = row[2].split(":")[0]
			ME_start = str(int(row[2].split(":")[1].split("-")[0]) -1)
			ME_end = row[2].split(":")[1].split("-")[1]
			ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end

			blencowe_ME.add(ME_ID)


	for row in csv.reader(open(ponting), delimiter = '\t'):

		if int(row[4])<=25:

			ME_chrom = row[0]
			ME_start = row[1]
			ME_end = row[2]
			ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end

			ponting_ME.add(ME_ID)

	mixture_cov = {}
	adipose_cov = {}
	adrenal_cov = {}
	brain_cov = {}
	breast_cov = {}
	colon_cov = {}
	heart_cov = {}
	kidney_cov = {}
	liver_cov = {}
	lung_cov = {}
	lymph_node_cov = {}
	ovary_cov = {}
	prostate_cov = {}
	skeletal_muscle_cov = {}
	testes_cov = {}
	thyroid_cov = {}
	white_blood_cells_cov = {}
	HepG2_control_cov = {}
	HepG2_UPF2_cov = {}
	HELA_control_cov = {}
	HELA_UPF1_cov = {}
	
	ME_SJ_coverage_reader(mixture, mixture_cov)
	ME_SJ_coverage_reader(adipose, adipose_cov)
	ME_SJ_coverage_reader(adrenal, adrenal_cov)
	ME_SJ_coverage_reader(brain, brain_cov)
	ME_SJ_coverage_reader(breast, breast_cov)
	ME_SJ_coverage_reader(colon, colon_cov)
	ME_SJ_coverage_reader(heart, heart_cov)
	ME_SJ_coverage_reader(kidney, kidney_cov)
	ME_SJ_coverage_reader(liver, liver_cov)
	ME_SJ_coverage_reader(lung, lung_cov)
	ME_SJ_coverage_reader(lymph_node, lymph_node_cov)
	ME_SJ_coverage_reader(ovary, ovary_cov)
	ME_SJ_coverage_reader(prostate, prostate_cov)
	ME_SJ_coverage_reader(skeletal_muscle, skeletal_muscle_cov)
	ME_SJ_coverage_reader(testes, testes_cov)
	ME_SJ_coverage_reader(thyroid, thyroid_cov)
	ME_SJ_coverage_reader(white_blood_cells, white_blood_cells_cov)
	ME_SJ_coverage_reader(HepG2_control, HepG2_control_cov)
	ME_SJ_coverage_reader(HepG2_UPF2, HepG2_UPF2_cov)
	ME_SJ_coverage_reader(HELA_control, HELA_control_cov)
	ME_SJ_coverage_reader(HELA_UPF1, HELA_UPF1_cov)

	for row in csv.reader(open(ME_centric_filter3), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_ME, total_ME, true_ME, score, is_annotated = row
		ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = true_ME.split("|")

		ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")
		ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end

		GENCODE = (is_annotated =="True")
		Blencowe = (ME_ID in blencowe_ME)
		Ponting = (ME_ID in ponting_ME)		


		is_alternative_5, is_alternative_3 = ME_alt[ME]

		# mixture_cov[ME]
		# adipose_cov[ME]
		# adrenal_cov[ME]
		# brain_cov[ME]
		# breast_cov[ME]
		# colon_cov[ME]
		# heart_cov[ME]
		# kidney_cov[ME]
		# liver_cov[ME]
		# lung_cov[ME]
		# lymph_node_cov[ME]
		# ovary_cov[ME]
		# prostate_cov[ME]
		# skeletal_muscle_cov[ME]
		# testes_cov[ME]
		# thyroid_cov[ME]
		# white_blood_cells_cov[ME]
		# HepG2_control_cov[ME]
		# HepG2_UPF2_cov[ME]
		# HELA_control_cov[ME]
		# HELA_UPF1_cov[ME]

		ME_cov_sum = 0
		SJ_cov_sum = 0
		
		for cov in [mixture_cov[ME], adipose_cov[ME], adrenal_cov[ME], brain_cov[ME], breast_cov[ME], colon_cov[ME], heart_cov[ME], kidney_cov[ME], liver_cov[ME], lung_cov[ME], lymph_node_cov[ME], ovary_cov[ME], prostate_cov[ME], skeletal_muscle_cov[ME], testes_cov[ME], thyroid_cov[ME], white_blood_cells_cov[ME], HepG2_control_cov[ME], HepG2_UPF2_cov[ME], HELA_control_cov[ME], HELA_UPF1_cov[ME]]:

			ME_cov, SJ_cov = cov.split("/")
			ME_cov = int(ME_cov)
			SJ_cov = sum(map(int, SJ_cov.split("|")))

			ME_cov_sum += ME_cov
			SJ_cov_sum += SJ_cov


		print ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, is_alternative_5, is_alternative_3, ME_cov_sum, SJ_cov_sum,  mixture_cov[ME], adipose_cov[ME], adrenal_cov[ME], brain_cov[ME], breast_cov[ME], colon_cov[ME], heart_cov[ME], kidney_cov[ME], liver_cov[ME], lung_cov[ME], lymph_node_cov[ME], ovary_cov[ME], prostate_cov[ME], skeletal_muscle_cov[ME], testes_cov[ME], thyroid_cov[ME], white_blood_cells_cov[ME], HepG2_control_cov[ME], HepG2_UPF2_cov[ME], HELA_control_cov[ME], HELA_UPF1_cov[ME] 

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23], sys.argv[24])


#python ~/my_src/ME/Pipeline/Round_2/Final_table.py /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 /media/HD3/Resultados/Micro_exons/Literature/Blencowe_et_al.txt.csv /media/HD3/Resultados/Micro_exons/Literature/Ponting_et_al.txt Mixture.TOTAL.filter1.ME_SJ_coverage


#python ~/my_src/ME/Pipeline/Round_2/Final_table.py /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 /media/HD3/Resultados/Micro_exons/Literature/Blencowe_et_al.txt.csv /media/HD3/Resultados/Micro_exons/Literature/Ponting_et_al.txt ../IBM/Mixture/Mixture.TOTAL.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/adipose.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/adrenal.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/brain.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/breast.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/colon.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/heart.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/kidney.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/liver.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/lung.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/lymph_node.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/ovary.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/prostate.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/skeletal_muscle.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/testes.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/thyroid.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../IBM/Tissue/TOTAL/white_blood_cells.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../shUPF1/TOTAL/HepG2_control.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../shUPF1/TOTAL/HepG2_shUPF2.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../shUPF1/TOTAL/sicontrol.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage ../shUPF1/TOTAL/siUPF1.TOTAL.sam.pre_processed.filter1.ME_SJ_coverage