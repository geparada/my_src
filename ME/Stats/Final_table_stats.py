import sys
import csv
from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes
from collections import defaultdict

def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][0], color='red')
    setp(bp['caps'][0], color='red')
    setp(bp['caps'][1], color='red')
    setp(bp['whiskers'][0], color='red')
    setp(bp['whiskers'][1], color='red')
    setp(bp['fliers'][0], color='red')
    setp(bp['fliers'][1], color='red')
    setp(bp['medians'][0], color='red')


def main(Final_ME_TABLE):

	######## BOXPLOT ############


	# novel_cov_len = defaultdict(list)
	# annotated_cov_len = defaultdict(list)


	# for row in csv.reader(open(Final_ME_TABLE), delimiter = ' '):


	# 	ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, ME_cov_sum, SJ_cov_sum,  mixture_cov[ME], adipose_cov[ME], adrenal_cov[ME], brain_cov[ME], breast_cov[ME], colon_cov[ME], heart_cov[ME], kidney_cov[ME], liver_cov[ME], lung_cov[ME], lymph_node_cov[ME], ovary_cov[ME], prostate_cov[ME], skeletal_muscle_cov[ME], testes_cov[ME], thyroid_cov[ME], white_blood_cells_cov[ME], HepG2_control_cov[ME], HepG2_UPF2_cov[ME], HELA_control_cov[ME], HELA_UPF1_cov[ME]  = row

	# 	len_micro_exon_seq_found = int(len_micro_exon_seq_found)


	# 	Novel = False

	# 	if GENCODE=="False" and Blencowe =="False" and Ponting=="False":

	# 		Novel = True

	# 	ME_cov_sum = 0
	# 	SJ_cov_sum = 0
		
	# 	for cov in [mixture_cov, adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov, HepG2_control_cov, HepG2_UPF2_cov, HELA_control_cov, HELA_UPF1_cov]:

	# 		ME_cov, SJ_cov = cov.split("/")
	# 		ME_cov = int(ME_cov)
	# 		SJ_cov = int(SJ_cov)

	# 		ME_cov_sum += ME_cov
	# 		SJ_cov_sum += SJ_cov_sum


	# 	if Novel:
	# 		novel_cov_len[len_micro_exon_seq_found].append(ME_cov_sum)

	# 	else:
	# 		annotated_cov_len[len_micro_exon_seq_found].append(ME_cov_sum)




	# fig = figure()
	# ax = axes()
	# hold(True)

	# C = 0

	# xticks_pos = []

	# for i in range(1,26):

	# 	data = [novel_cov_len[i], annotated_cov_len[i]]
	# 	bp = boxplot(data, positions = [C + i, C + i + 1], widths=0.6)
	# 	setBoxColors(bp)

	# 	xticks_pos.append(C + i + 0.5) 

	# 	C += 3

	# #xlim(3,78)
	# ylim(0,100)

	# ax.set_xticklabels(map(str,range(1,26)))
	# ax.set_xticks(xticks_pos)

	# # draw temporary red and blue lines and use them to create a legend
	# hB, = plot([1,1],'b-')
	# hR, = plot([1,1],'r-')
	# legend((hR, hB),('Novel', 'Annotated'))
	# hB.set_visible(False)
	# hR.set_visible(False)

	# savefig('Novel_Annotated_coverage.png')
	#show()

	####### Tisue specific ###########

	

	for row in csv.reader(open(Final_ME_TABLE), delimiter = ' '):


		ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, ME_cov_sum, SJ_cov_sum,  mixture_cov[ME], adipose_cov[ME], adrenal_cov[ME], brain_cov[ME], breast_cov[ME], colon_cov[ME], heart_cov[ME], kidney_cov[ME], liver_cov[ME], lung_cov[ME], lymph_node_cov[ME], ovary_cov[ME], prostate_cov[ME], skeletal_muscle_cov[ME], testes_cov[ME], thyroid_cov[ME], white_blood_cells_cov[ME], HepG2_control_cov[ME], HepG2_UPF2_cov[ME], HELA_control_cov[ME], HELA_UPF1_cov[ME] = row

		len_micro_exon_seq_found = int(len_micro_exon_seq_found)

		ME_covs = []
		SJ_covs = []

		for cov in [mixture_cov, adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov]:

			ME_cov, SJ_cov = cov.split("/")
			ME_cov = int(ME_cov)
			SJ_cov = int(SJ_cov)

			ME_covs.append(ME_cov)
			SJ_covs.append(SJ_cov)



		if Novel:
			novel_cov_len[len_micro_exon_seq_found].append(ME_cov_sum)

		else:
			annotated_cov_len[len_micro_exon_seq_found].append(ME_cov_sum)



if __name__ == '__main__':
	main(sys.argv[1])
