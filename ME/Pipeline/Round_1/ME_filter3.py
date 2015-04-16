import sys
import csv
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np


def xfrange(start, stop, step):
    while start < stop:
        yield start
        start += step



def main(ME_centric, simulation_scores, gencode_bed):

	simulation_U2 = []
	simulation_vertebrates = []
	simulation_primates = []

	gencode_exons = set([])

	for row in csv.reader(open(gencode_bed), delimiter = '\t'):

		csv.field_size_limit(1000000000)

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]


		for q1, b in zip(qstarts[1:-1], blocksizes[1:-1]):
			estart = start + q1
			eend = start + q1 + b

			exon = "_".join(map(str, [chr, strand, estart, eend]))

			gencode_exons.add(exon)




	for row in csv.reader(open(simulation_scores), delimiter = ' '):

		chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates = row

		simulation_U2.append(float(U2_score))
		simulation_vertebrates.append(float(mean_conservation_vertebrates))
		simulation_primates.append(float(mean_conservation_primates))


	for row in csv.reader(open(ME_centric), delimiter = '\t'):

		annotated = False

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, P_ME, total_ME = row

		total_number_of_micro_exons_matches = float(total_number_of_micro_exons_matches)

		true_ME = []
		len_micro_exon_seq_found = int(len_micro_exon_seq_found)

		if total_number_of_micro_exons_matches <= 5 and len_micro_exon_seq_found <= 25:

			for m in total_ME.split(","):

				if m!="":

					ME, U2_score, mean_conservation_vertebrates, mean_conservation_primates = m.split("|")

					U2_FPR = 1 - stats.percentileofscore(simulation_U2, float(U2_score))/100
					vertebrates_FPR = 1 - stats.percentileofscore(simulation_vertebrates, float(mean_conservation_vertebrates))/100
					primates_FPR = 1 - stats.percentileofscore(simulation_primates, float(mean_conservation_primates))/100

					#conservation_FPR = min([vertebrates_FPR, primates_FPR])
					conservation_FPR = vertebrates_FPR

					total_score = U2_FPR * conservation_FPR * float(total_number_of_micro_exons_matches)

					if total_score <= 0.05 and U2_FPR <= 0.2:
						if len_micro_exon_seq_found >= 7:
							true_ME.append((m, total_score))
						elif total_score <= 0.005 and vertebrates_FPR <= 0.05 and U2_FPR <= 0.2:
							true_ME.append((m, total_score))


		if len(true_ME)==1: #and int(len_micro_exon_seq_found)>25:

			

			if true_ME[0][0].split("|")[0] in gencode_exons:

				annotated = True
				#print true_ME[0][0].split("|")[0]

			print "\t".join(row) + "\t" + "\t".join(map(str, true_ME[0])) + "\t" + str(annotated) 

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])	


#python ~/my_src/ME/Pipeline/ME_filter3.py TOTAL.filter1.ME_centric.filter2 ../../simulation/ROC/ME_sim.scores ~/db/transcriptome/hg19/Gene_models/gencode/v17/gencode.v17.annotation.bed12