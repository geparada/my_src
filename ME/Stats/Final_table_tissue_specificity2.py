import sys
import csv
from collections import defaultdict
import numpy as np
import math


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

	out = []

	for row in csv.reader(open(Final_ME_TABLE), delimiter = ' '):


		ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, ME_cov_sum, SJ_cov_sum,  mixture_cov, adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov, HepG2_control_cov, HepG2_UPF2_cov, HELA_control_cov, HELA_UPF1_cov = row

		len_micro_exon_seq_found = int(len_micro_exon_seq_found)

		ME_covs = []
		SJ_covs = []

		all_tissues_phi = []

		tissues_phi = []
		tissues_number = []

		C = 0

		for cov in [adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov]:

			ME_cov, SJ_cov = cov.split("/")
			ME_cov = int(ME_cov)
			SJ_cov = int(SJ_cov)
			Total_cov = ME_cov + SJ_cov


			if Total_cov >= 10:

				#phi, b_low, b_up = calcBin(ME_cov , Total_cov)
				# ME_covs.append(ME_cov)
				# SJ_covs.append(SJ_cov)				

				phi = float(ME_cov)/float(Total_cov)

				tissues_phi.append(phi)
				tissues_number.append(C)

			all_tissues_phi.append(phi)
			C += 1

		if np.std(tissues_phi)!=0:

			tissues_phi_norm = []

			for i in tissues_phi:
				phi_norm = float(i - np.mean(tissues_phi)) / float(np.std(tissues_phi))
				tissues_phi_norm.append(phi_norm)

			phi_norm_tissue_number = []

			for p, t in zip(tissues_phi_norm, tissues_number):
				phi_norm_tissue_number.append([p, t])

			phi_norm_tissue_number.sort(key=lambda x:x[0], reverse=True)

			Us = []

			for i in range(1,len(phi_norm_tissue_number)):

				outilers_phi = []
				outilers_tissue_number = []
				non_outilers_phi = []
				non_outilers_tissue_number = []	

				for p, t in phi_norm_tissue_number[:i]:

					outilers_phi.append(p)
					outilers_tissue_number.append(t)

				for p, t in phi_norm_tissue_number[i:]:

					non_outilers_phi.append(p)
					non_outilers_tissue_number.append(t)

				sigma = np.std(non_outilers_phi)

				# if sigma==0:
				# 	sigma = 0.0001
				n = len(non_outilers_phi)
				s = len(outilers_phi)

				U = float(np.mean(outilers_phi) - np.mean(non_outilers_phi)) / float(np.std(non_outilers_phi)  + 1)
				#U = n * np.log(sigma) + float(np.sqrt(2)*s*np.log(np.math.factorial(n)))/float(n)

				Us.append((U, outilers_tissue_number, non_outilers_tissue_number))


			if len(Us)!=0:
			

				Umax, outilers_tissue_number, non_tissue_number = max(Us, key=lambda x: x[0])

				Umax_phi = []
				Umax_tissues = []
				Umax_non_phi = []
				Umax_non_tissues = []

				for n in outilers_tissue_number:

					#print len(tissues_phi), n
					#print tissues_phi[n]
					Umax_phi.append(all_tissues_phi[n])
					Umax_tissues.append(tissues[n])

				for n in non_tissue_number:


					Umax_non_phi.append(all_tissues_phi[n])
					Umax_non_tissues.append(tissues[n])


				TE_score = float((np.mean(Umax_phi)-np.mean(Umax_non_phi)) * len(Umax_non_phi)) / float(len(Umax_phi) + len(Umax_non_phi))




				out.append([ ME, ME_cov_sum, SJ_cov_sum, len_micro_exon_seq_found, micro_exon_seq_found, GENCODE, Blencowe, Ponting, Umax, TE_score, ",".join(Umax_tissues), Umax_phi, Umax_non_phi ])



	out.sort(key=lambda x:x[9], reverse=True)

	for i in out:

		print "\t".join(map(str, i))





if __name__ == '__main__':
	main(sys.argv[1])







