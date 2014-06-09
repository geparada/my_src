import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import scipy.stats


IE_kmers = defaultdict(set)

total_non_can = defaultdict(int)
non_can_EI = defaultdict(int)
total_can = defaultdict(int)
can_EI = defaultdict(int)

def index (c, total):
	try:
		return (float(c))/float(total)
	except ZeroDivisionError:
		return 0

def main (non_can_fasta, can_fasta, EIE_cluster, direction):

	#direction = 53 (para intron5, exon3) or 35 (para exon5, intron3)

	for row in csv.reader(open(EIE_cluster), delimiter = '\t'):
		kmer, info = row[:2]   #Es necesario ya que el archivo tiene unos espacios demas a la derecha en algunas filas
		IE, SRE, sub_SRE = info.split("/")
		IE_kmers[kmer].add(IE)

	for record in SeqIO.parse(non_can_fasta, "fasta"):
		kmer_len = 6
		seq = str(record.seq)

		if direction == '35':
			seq = str(record.seq)[::-1]			
			
		for i in range(len(seq)- kmer_len):
			kmer = str(seq[i:i+kmer_len])
			IE_types = []
				
			if IE_kmers.has_key(kmer):
				IE_types = list(IE_kmers[kmer])
				for t in IE_types:
					for pos in range(i, i+kmer_len):
						event_ID = t + "_" + str(pos)   # el event_ID indica que EIE es y en que base se encontro.
						non_can_EI[event_ID] += 1

			for pos in range(i, i+kmer_len):
				total_non_can[pos] += 1						

	for record in SeqIO.parse(can_fasta, "fasta"):
		kmer_len = 6
		seq = str(record.seq)

		if direction == '35':
			seq = str(record.seq)[::-1]			
			
		for i in range(len(seq)- kmer_len):
			kmer = str(seq[i:i+kmer_len])
			IE_types = []
				
			if IE_kmers.has_key(kmer):
				IE_types = list(IE_kmers[kmer])
				for t in IE_types:
					for pos in range(i, i+kmer_len):
						event_ID = t + "_" + str(pos)   # el event_ID indica que EIE es y en que base se encontro.
						can_EI[event_ID] += 1

			for pos in range(i, i+kmer_len):
				total_can[pos] += 1
									
	print "nt", "index_non_can_EIE", "index_can_EIE", "index_non_can_IIE", "index_can_IIE"
	
	c = 0
	count_index_non_can_EIE = 0
	count_index_can_EIE = 0
	count_index_non_can_IIE = 0
	count_index_can_IIE = 0
	
	interval = 10
	L = 100
	


	for i in range(L):
		count_non_can_total = total_non_can[i]
		count_can_total = total_can[i]	
		
		index_non_can_EIE = index(non_can_EI["EIE_" + str(i)], count_non_can_total)
		index_can_EIE = index(can_EI["EIE_" + str(i)], count_can_total)
		
		index_non_can_IIE = index(non_can_EI["IIE_" + str(i)], count_non_can_total)
		index_can_IIE = index(can_EI["IIE_" + str(i)], count_can_total)
		
		c += 1
		
		if c<interval:
			count_index_non_can_EIE += index_non_can_EIE
			count_index_can_EIE += index_can_EIE
			count_index_non_can_IIE += index_non_can_IIE
			count_index_can_IIE += index_can_IIE
		
		if c==interval:
			mean_index_non_can_EIE = index(count_index_non_can_EIE, interval )
			mean_index_can_EIE = index(count_index_can_EIE, interval )
			mean_index_non_can_IIE = index(count_index_non_can_IIE, interval )
			mean_index_can_IIE = index(count_index_can_IIE, interval )
			
			print i+1, mean_index_non_can_EIE, mean_index_can_EIE, mean_index_non_can_IIE, mean_index_can_IIE
						
			c=0
			count_index_non_can_EIE = 0
			count_index_can_EIE = 0
			count_index_non_can_IIE = 0
			count_index_can_IIE = 0




#python ~/my_src/Motifs/SJ_graph.py ../no_canonicos5_intron5 ../canonicos_intron5 TOTAL_clusters_5_IE 53

#python ~/my_src/Motifs/SJ_graph.py ../no_canonicos5_exon3 ../canonicos_exon3 TOTAL_clusters_5_IE 53


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])	


