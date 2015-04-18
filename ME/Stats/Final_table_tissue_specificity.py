import sys
import csv
from collections import defaultdict
import numpy as np


def U_combinations(S0, S1, c, R):

	if c==1:
		for a in S0:
			R.append(a)

	if c < len(S0) - 1:

		S2 = []

		for a in S1:

			for b in S0[max(map(int,a.split("-")))+1:]:				
				S2.append(a + "-" + b)
				R.append(a + "-" + b)

		U_combinations(S0, S2, c+1, R)

	else:
		return


def binP(N, p, x1, x2):
    p = float(p)
    q = p/(1-p)
    k = 0.0
    v = 1.0
    s = 0.0
    tot = 0.0

    while(k<=N):
            tot += v
            if(k >= x1 and k <= x2):
                    s += v
            if(tot > 10**30):
                    s = s/10**30
                    tot = tot/10**30
                    v = v/10**30
            k += 1
            v = v*q*(N+1-k)/k
    return s/tot

def calcBin(vx, vN, vCL = 95):
	'''
	Calculate the exact confidence interval for a binomial proportion

	Usage:
	>>> calcBin(13,100)    
	(0.07107391357421874, 0.21204372406005856)
	>>> calcBin(4,7)   
	(0.18405151367187494, 0.9010086059570312)
	'''


 
	vx = float(vx)
	vN = float(vN)
	#Set the confidence bounds
	vTU = (100 - float(vCL))/2
	vTL = vTU

	vP = vx/vN
	if(vx==0):
		dl = 0.0
	else:
		v = vP/2
		vsL = 0
		vsH = vP
		p = vTL/100

		while((vsH-vsL) > 10**-5):
			if(binP(vN, v, vx, vN) > p):
				vsH = v
				v = (vsL+v)/2
			else:
				vsL = v
				v = (v+vsH)/2
		dl = v

	if(vx==vN):
		ul = 1.0
	else:
		v = (1+vP)/2
		vsL =vP
		vsH = 1
		p = vTU/100
		while((vsH-vsL) > 10**-5):
			if(binP(vN, v, 0, vx) < p):
				vsH = v
				v = (vsL+v)/2
			else:
				vsL = v
				v = (v+vsH)/2
		ul = v

	return vP, vP - dl, ul - vP          
	#print '%s\t%s\t%s' % (vP, vP - dl, ul - vP)


def main(Final_ME_TABLE):

	tissues = {0:"adipose", 1:"adrenal", 2:"brain", 3:"breast", 4:"colon", 5:"heart", 6:"kidney", 7:"liver", 8:"lung", 9:"lymph_node", 10:"ovary", 11:"prostate", 12:"skeletal_muscle", 13:"testes", 14:"thyroid", 15:"white_blood_cells"}


	for row in csv.reader(open(Final_ME_TABLE), delimiter = ' '):


		ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, ME_cov_sum, SJ_cov_sum,  mixture_cov, adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov, HepG2_control_cov, HepG2_UPF2_cov, HELA_control_cov, HELA_UPF1_cov = row

		len_micro_exon_seq_found = int(len_micro_exon_seq_found)

		ME_covs = []
		SJ_covs = []

		tissues_phi = []
		tissues_number = []

		C = 0

		for cov in [adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov]:

			ME_cov, SJ_cov = cov.split("/")
			ME_cov = int(ME_cov)
			SJ_cov = int(SJ_cov)
			Total_cov = ME_cov + SJ_cov


			if Total_cov != 0:

				#phi, b_low, b_up = calcBin(ME_cov , Total_cov)
				# ME_covs.append(ME_cov)
				# SJ_covs.append(SJ_cov)				

				phi = float(ME_cov)/float(Total_cov)

				tissues_phi.append(phi)
				tissues_number.append(C)

			C += 1

		if np.std(tissues_phi)!=0:

			tissues_phi_norm = []

			for i in tissues_phi:
				phi_norm = float(i - np.mean(tissues_phi)) / float(np.std(tissues_phi))
				tissues_phi_norm.append(phi_norm)


			non_outilers_com = []

			list_U_com = map(str, range(len(tissues_phi)))

			U_combinations(list_U_com, list_U_com, 1, non_outilers_com)

			Us = []


			for com in non_outilers_com:

				non_outilers_phi = []
				non_outliers = map(int, com.split("-"))
				outliers = set(range(len(tissues_phi_norm))) - set(non_outliers)

				for i in non_outliers:
					non_outilers_phi.append(tissues_phi_norm[i])

				sigma = np.std(non_outilers_phi)

				n = len(non_outliers)
				s = len(tissues_phi_norm) - n

				#U = 10**10000000000

				if sigma != 0:

					U = n * np.log(sigma) + float(np.sqrt(2)*s*np.log(np.math.factorial(n)))/float(n)

					Us.append([U, outliers, non_outliers])

			if len(Us)!=0:
			
				Umin, Umin_outliers, Umin_non_outliers  = min(Us, key=lambda x: x[0])

				Umin_phi = []
				Umin_tissues = []
				Umin_non_phi = []
				Umin_non_tissues = []

				for n in Umin_outliers:

					Umin_phi.append(tissues_phi[n])
					Umin_tissues.append(tissues[tissues_number[n]])

				for n in Umin_non_outliers:

					Umin_non_phi.append(tissues_phi[n])
					Umin_non_tissues.append(tissues[tissues_number[n]])


				if np.mean(Umin_phi) >  np.mean(Umin_non_phi):

					print ME, ME_cov_sum, SJ_cov_sum, len_micro_exon_seq_found, micro_exon_seq_found, GENCODE, Blencowe, Ponting, Umin, np.mean(Umin_phi)-np.mean(Umin_non_phi), Umin_tissues, Umin_phi, Umin_non_phi

				else:
					print ME, ME_cov_sum, SJ_cov_sum, len_micro_exon_seq_found, micro_exon_seq_found, GENCODE, Blencowe, Ponting, Umin, np.mean(Umin_non_phi)-np.mean(Umin_phi), Umin_non_tissues, Umin_non_phi, Umin_phi








if __name__ == '__main__':
	main(sys.argv[1])







