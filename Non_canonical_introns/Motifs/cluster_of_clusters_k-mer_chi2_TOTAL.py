import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import scipy.stats

non_can_kmers = defaultdict(int)
can_kmers = defaultdict(int)

clusters_non_can = defaultdict(int)
clusters_can = defaultdict(int)

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
				

def main (non_can_fasta, can_fasta, kmer_clusters):
	""" Realiza una tabla de frecuencias de los k-meks precentes """

	non_can_seq_count = 0
	can_seq_count = 0

	non_can_seq_len_count = 0
	can_seq_len_count = 0


	results = []
	
	kmer_lens = set([])         #Conjunto de todos los tamanos de cluster que hay
	cluster_kmer_lens = {}

	high_cluster = []

	for row in csv.reader(open(kmer_clusters), delimiter = '\t'):
		kmer = row[0]
		cluster_id = row[1].split("/")[0]
		kmer_lens.add(len(kmer))
		cluster_kmer_lens[cluster_id] = len(kmer)


	for record in SeqIO.parse(non_can_fasta, "fasta"):

		non_can_seq_count += 1
		non_can_seq_len_count += len(record)

		for kmer_len in kmer_lens:     #Para cada tamano de clusters se extraen los k-mer
			for i in range(len(record.seq)- kmer_len):
				non_can_kmers[str(record.seq[i:i+kmer_len])] += 1

	for record in SeqIO.parse(can_fasta, "fasta"):
		
		can_seq_count += 1
		can_seq_len_count += len(record)		
		
		for kmer_len in kmer_lens:             
			for i in range(len(record.seq)- kmer_len):
				can_kmers[str(record.seq[i:i+kmer_len])] += 1

	for row in csv.reader(open(kmer_clusters), delimiter = '\t'):
		kmer = row[0]
		cluster_id = row[1].split("/")[0]
		
		non_can_count = non_can_kmers[kmer]
		can_count = can_kmers[kmer]
		clusters_non_can[cluster_id] += non_can_count
		clusters_can[cluster_id] += can_count
		cluster_key.add(cluster_id)

	for k in cluster_key:
		kmer_len = cluster_kmer_lens[k]
		
		non_can_nkmers = non_can_seq_len_count  -  (kmer_len * non_can_seq_count)
		can_nkmers = can_seq_len_count  -  (kmer_len * can_seq_count)
		
		non_can_count = clusters_non_can[k]
		can_count = clusters_can[k]		
		expected_non_can_count = (float(non_can_nkmers) / float(can_nkmers)) * can_count		
		chi2 = chi_square(non_can_count, expected_non_can_count)
		
		results.append((k, non_can_count, expected_non_can_count, chi2))
	
	GL = len(results) -1 
	
	for r in sorted(results, key=lambda a: a[3], reverse=True):
		cluster = r[0]
		non_can = r[1]
		expected_non_can = r[2]
		chi2 = r[3]
		p_value = scipy.stats.chi2.sf(chi2, GL)
		
		print cluster, non_can, expected_non_can, chi2, p_value
	



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])	



