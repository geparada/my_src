import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0

def main (non_can_fasta, can_fasta, kmer_clusters):
	""" Cuenta el numero de secuencias que hace match con los clusters """
	
	clusters = defaultdict(list)
	non_can_clusters_count = defaultdict(int)
	can_clusters_count = defaultdict(int)
	cluster_key = set([])
	total_non_can = 0
	total_can = 0
	
	
	for row in csv.reader(open(kmer_clusters), delimiter = '\t'):
		kmer = row[0]
		cluster_id = row[1]
		clusters[cluster_id].append(kmer)
		cluster_key.add(cluster_id)



	for record in SeqIO.parse(non_can_fasta, "fasta"):
		
		total_non_can += 1
		
		for c in clusters.items():
			cluster_id = c[0]
			kmer = c[1]
			match = False
			
			for k in kmer:
				if k in record.seq:
					match = True
				if match == True:
					break
			
			if match == True:
				non_can_clusters_count[cluster_id] += 1



	for record in SeqIO.parse(can_fasta, "fasta"):
		
		total_can += 1
		
		for c in clusters.items():
			cluster_id = c[0]
			kmer = c[1]
			match = False
			
			for k in kmer:
				if k in record.seq:
					match = True
				if match == True:
					break
			
			if match == True:
				can_clusters_count[cluster_id] += 1


	for k in cluster_key:
		non_can = percent(non_can_clusters_count[k], total_non_can)
		can = percent(can_clusters_count[k], total_can)
		
		
		print k, non_can, can, (non_can + can) / float(2)	
		


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])	
