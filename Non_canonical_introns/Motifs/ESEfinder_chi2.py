
import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import scipy.stats

non_can_hexamers = defaultdict(int)
can_hexamers = defaultdict(int)
ESE_clusters_non_can = defaultdict(int)
ESE_clusters_can = defaultdict(int)
cluster_key = set([])

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0
		
def chi_square (o, e):
	""" Calculate the chi-square """
	try:
		return ((o -e)**2)/e
	except ZeroDivisionError:
		return 0		

def k_mer_analizer (seq, matrix, k_mer_dict, k_mer_count):
	kmer_len = matrix["len"]
	threshold = matrix["threshold"]
	k_mers = []
	
	for i in range(len(seq)- kmer_len):
		k_mer = str(seq[i:i+kmer_len]).upper()
		
		for n in k_mer
		
	
	
	k_mer_dict[str(seq[i:i+kmer_len])]
	
	
	""" extrae y analiza k-mers de los fasta entregados"""
					

def main (non_can_fasta, can_fasta):
	""" Analiza si hay sobrerrepesentacion de secuencias que cumplan con las matrices """

	SRSF1A = {"len":7, "threshold":1.956, "A":[-1.14, 0.62, -1.58, 1.32, -1.58, -1.58, 0.62], "C":[1.37, -1.1, 0.73, 0.33, 0.94, -1.58, -1.58], "G":[-0.21, 0.17, 0.48, -1.58, 0.33, 0.99, -0.11], "T":[-1.58, -0.5, -1.58, -1.13, -1.58, -1.13, 0.27]}
	SRSF1B = {"len":7, "threshold":1.867, "A":[-1.58, 0.15, -0.97, 0.74, -1.19, -0.75, 0.43], "C":[1.55, -0.53, 0.79, 0.33, 0.72, -0.62, -0.99], "G":[-1.35, 0.44, 0.41, -0.98, 0.51, 1.03, 0.00], "T":[-1.55, -0.28, -1.28, -0.92, -1.09, -0.52, 0.20]}
	SRSF2 =	{"len":8, "threshold":2.383, "A":[-0.88, 0.09, -0.06, -1.58, 0.09, -0.41, -0.06, 0.23], "C":[1.16, -1.58, 0.95, 1.11, 0.56, 0.86, 0.32, -1.58], "G":[0.87, 0.45, -1.36, -1.58, -0.33, -0.05, -1.36, 0.68], "T":[-1.18, -0.2, 0.38, 0.88, -0.2, -0.86, 0.96, -1.58]}
	SRSF5 = {"len":7, "threshold":2.670, "A":[-0.13, -1.58, 1.28, -0.33, 0.97, -0.13, -1.58], "C":[0.56, 0.68, -1.12, 1.24, -0.77, 0.13, -0.05], "G":[-1.58, -0.14, -1.33, -0.48, -1.58, 0.44, 0.8], "T":[0.92, 0.37, 0.23, -1.14, 0.72, -1.58, -1.58]}
	SRSF6 = {"len":6, "threshold":2.676, "A":[-0.66, 0.11, -0.66, 0.11, -1.58, 0.61], "C":[0.39, -1.58, 1.48, -1.58, -1.58, 0.98], "G":[-1.58, 0.72, -1.58, 0.72, 0.21, -0.79], "T":[1.22, -1.58, -0.07, -1.58, 1.02, -1.58]}

	reader = csv.reader(open(k_mek_clusters), delimiter = '\t')

	non_can_nhexamers = 0
	can_nhexamers = 0
	
	results = []

	for record in SeqIO.parse(non_can_fasta, "fasta"):
#		kmer_len = 6
		
		for i in range(len(record.seq)- kmer_len):
			non_can_hexamers[str(record.seq[i:i+kmer_len])] += 1
			non_can_nhexamers += 1

	for record in SeqIO.parse(can_fasta, "fasta"):
#		kmer_len = 6
		
		for i in range(len(record.seq)- kmer_len):
			can_hexamers[str(record.seq[i:i+kmer_len])] += 1
			can_nhexamers += 1		

	for row in reader:
		hexamer = row[0]
		cluster_id = row[1]
		
		non_can_count = non_can_hexamers[hexamer]
		can_count = can_hexamers[hexamer]
		ESE_clusters_non_can[cluster_id] += non_can_count
		ESE_clusters_can[cluster_id] += can_count
		cluster_key.add(cluster_id)		


	for k in cluster_key:
		non_can_count = ESE_clusters_non_can[k]
		can_count = ESE_clusters_can[k]		
		expected_non_can_count = (float(non_can_nhexamers) / float(can_nhexamers)) * can_count		
		chi2 = chi_square(non_can_count, expected_non_can_count)
		results.append((k, non_can_count, expected_non_can_count, chi2))
	
	
	GL = len(results) -1 
	
	for r in sorted(results, key=lambda a: a[3], reverse=True):
		hexamer = r[0]
		non_can = r[1]
		expected_non_can = r[2]
		chi2 = r[3]		
		p_value = scipy.stats.chi2.sf(chi2, GL)
		
		print hexamer, non_can, expected_non_can, chi2, p_value
	



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))	
