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
				

def main (non_can_fasta, can_fasta, ESE_hexamers):
	""" Realiza una tabla de frecuencias de los k-meks precentes """

	reader = csv.reader(open(ESE_hexamers), delimiter = '\t')

	non_can_nhexamers = 0
	can_nhexamers = 0
	
	results = []

	for record in SeqIO.parse(non_can_fasta, "fasta"):
		kmer_len = 6
		
		for i in range(len(record.seq)- kmer_len):
			non_can_hexamers[str(record.seq[i:i+kmer_len])] += 1
			non_can_nhexamers += 1

	for record in SeqIO.parse(can_fasta, "fasta"):
		kmer_len = 6
		
		for i in range(len(record.seq)- kmer_len):
			can_hexamers[str(record.seq[i:i+kmer_len])] += 1
			can_nhexamers += 1		

	for row in reader:
		hexamer = row[0]
		cluster_id = row[1]
		
		non_can_count = non_can_hexamers[hexamer]
		can_count = can_hexamers[hexamer]
#		expected_non_can = (float(non_can_nhexamers) / float(can_nhexamers)) * can_count		
#		chi2 = chi_square(non_can_count, expected_non_can)
		ESE_clusters_non_can[cluster_id] += non_can_count
		ESE_clusters_can[cluster_id] += can_count
		cluster_key.add(cluster_id)		
		
#		results.append((hexamer, non_can, expected_non_can, chi2))


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
	main(sys.argv[1], sys.argv[2], sys.argv[3] )	



