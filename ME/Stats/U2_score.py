import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def PWM_to_dict(file):
	reader = csv.reader(open(file), delimiter = '\t')
	header = reader.next()
	header_dict = {}
	col = 0
	
	matrix = {}
	
	for name in header:
		header_dict[name] = col
		col += 1
	
	A_frec = []
	C_frec = []
	G_frec = []
	T_frec = []
	
	for row in reader:
		A = row[header_dict["A"]]
		C = row[header_dict["C"]]
		G = row[header_dict["G"]]
		T = row[header_dict["T"]]
		
		A_frec.append(float(A))
		C_frec.append(float(C))
		G_frec.append(float(G))
		T_frec.append(float(T))
			
	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	matrix["N"] = min(A_frec , C_frec, G_frec, T_frec)
	
	return matrix


def main(ME_final, U2_GTAG_5_file, U2_GTAG_3_file):

	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])



	for row in csv.reader(open(ME_final), delimiter = '\t'):

		#sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_ME, total_ME, true_ME, score, is_annotated = row

		#ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, ME_cov_sum, SJ_cov_sum,  mixture_cov, adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov, HepG2_control_cov, HepG2_UPF2_cov, HELA_control_cov, HELA_UPF1_cov = row

		ME_chrom, ME_start,ME_end, ME, len_micro_exon_seq_found, ME_strand = row


		#ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")


		ME_start = int(ME_start)
		ME_end = int(ME_end)


		ME5 = str(Genome[ME_chrom][ME_start-14:ME_start+3]).upper()
		ME3 = str(Genome[ME_chrom][ME_end-3:ME_end+10]).upper()


		if ME_strand == "-":

			ME5 = str(Genome[ME_chrom][ME_end-3:ME_end+14].reverse_complement()).upper()
			ME3 = str(Genome[ME_chrom][ME_start-10:ME_start+3].reverse_complement()).upper()


		U2_score_5 = 0
		U2_score_3 = 0

		i = 0

		for N in ME5:
			U2_score_3 += U2_GTAG_3[N][i]
			i += 1

		i = 0

		for N in ME3:
			U2_score_5 += U2_GTAG_5[N][i]
			i += 1

		U2_score_5 = percent(U2_score_5, U2_GTAG_5_max_score)
		U2_score_3 = percent(U2_score_3, U2_GTAG_3_max_score)


		print ME, len_micro_exon_seq_found, U2_score_3, U2_score_5


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2], sys.argv[3], sys.argv[4])


# python ~/my_src/ME/Stats/U2_score.py ~/db/genome/hg19.fa Final_ME_table.5cov.3len_ME.5alt ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix