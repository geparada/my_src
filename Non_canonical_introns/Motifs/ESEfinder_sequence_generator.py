import sys
import csv
import math


def sequence_generator (L):
	""" Genera todas las posibles sequencias de un largo L """
	
	seqs = ["A", "C", "G", "T"]
	
	while L > 1:
		L -= 1
		new_seqs = []
		
		for a in seqs:
			for b in "ACGT":
				new_seqs.append(a + b)
		
		seqs = new_seqs
	
	return seqs			 

def PWM_to_seqs (kmer, matrix):
	""" Calcula el conjunto de secuencias probables que podrian superar el threshold para cada matriz """
	
	kmer_len = matrix["len"]
	threshold = matrix["threshold"]
	cluster_name =  matrix["name"]

	score = 0

	for i, n in zip(range(kmer_len), kmer):
		score += matrix[n][i]

	if score > threshold:
		print "%s\t%s" % (kmer, cluster_name)

	

def main ():
	""" Obtiene lista de secuencias que pertenecen a cada motivo del ESEfinder """

	
	SRSF1A = {"name":"SRSF1A", "len":7, "threshold":1.956, "A":[-1.14, 0.62, -1.58, 1.32, -1.58, -1.58, 0.62], "C":[1.37, -1.1, 0.73, 0.33, 0.94, -1.58, -1.58], "G":[-0.21, 0.17, 0.48, -1.58, 0.33, 0.99, -0.11], "T":[-1.58, -0.5, -1.58, -1.13, -1.58, -1.13, 0.27]}
	SRSF1B = {"name":"SRSF1B", "len":7, "threshold":1.867, "A":[-1.58, 0.15, -0.97, 0.74, -1.19, -0.75, 0.43], "C":[1.55, -0.53, 0.79, 0.33, 0.72, -0.62, -0.99], "G":[-1.35, 0.44, 0.41, -0.98, 0.51, 1.03, 0.00], "T":[-1.55, -0.28, -1.28, -0.92, -1.09, -0.52, 0.20]}
	SRSF2 =	{"name":"SRSF2", "len":8, "threshold":2.383, "A":[-0.88, 0.09, -0.06, -1.58, 0.09, -0.41, -0.06, 0.23], "C":[1.16, -1.58, 0.95, 1.11, 0.56, 0.86, 0.32, -1.58], "G":[0.87, 0.45, -1.36, -1.58, -0.33, -0.05, -1.36, 0.68], "T":[-1.18, -0.2, 0.38, 0.88, -0.2, -0.86, 0.96, -1.58]}
	SRSF5 = {"name":"SRSF5", "len":7, "threshold":2.670, "A":[-0.13, -1.58, 1.28, -0.33, 0.97, -0.13, -1.58], "C":[0.56, 0.68, -1.12, 1.24, -0.77, 0.13, -0.05], "G":[-1.58, -0.14, -1.33, -0.48, -1.58, 0.44, 0.8], "T":[0.92, 0.37, 0.23, -1.14, 0.72, -1.58, -1.58]}
	SRSF6 = {"name":"SRSF6", "len":6, "threshold":2.676, "A":[-0.66, 0.11, -0.66, 0.11, -1.58, 0.61], "C":[0.39, -1.58, 1.48, -1.58, -1.58, 0.98], "G":[-1.58, 0.72, -1.58, 0.72, 0.21, -0.79], "T":[1.22, -1.58, -0.07, -1.58, 1.02, -1.58]}

	hexamers = sequence_generator(6)
	heptamers = sequence_generator(7)	
	octamers = sequence_generator(8)	
	
	for k in hexamers:
		PWM_to_seqs(k, SRSF6)

	for k in heptamers:
		PWM_to_seqs(k, SRSF1A)
		PWM_to_seqs(k, SRSF1B)
		PWM_to_seqs(k, SRSF5)

	for k in octamers:
		PWM_to_seqs(k, SRSF2)
		
	

main()




#if __name__ == '__main__':
#	sequence_generator(int(sys.argv[1]))
