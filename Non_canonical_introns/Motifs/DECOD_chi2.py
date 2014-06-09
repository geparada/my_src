import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import scipy.stats


def chi_square (o, e):
	""" Calculate the chi-square """
	try:
		return ((o -e)**2)/e
	except ZeroDivisionError:
		return 0		

def main (DECOD_OUTPUT, non_can_fasta, can_fasta):
	""" Function doc """
	
	kmer_len = 8
	
	pos_neg = -1
	n_motif = 0
	pos_motifs = defaultdict(int)
	neg_motifs = defaultdict(int)
	keys = set([])
	
	non_can_seq_count = 0
	non_can_seq_len_count = 0
	can_seq_count = 0
	can_seq_len_count = 0
	
	results = []
	
	
	for record in SeqIO.parse(non_can_fasta, "fasta"):

		non_can_seq_count += 1
		non_can_seq_len_count += len(record)

	for record in SeqIO.parse(can_fasta, "fasta"):
		
		can_seq_count += 1
		can_seq_len_count += len(record)		
		
	
	for row in csv.reader(open(DECOD_OUTPUT), delimiter = ' '):
		
		try:	
			
			if row[0][:6] == ">Motif":
				n_motif += 1
				keys.add("Motif"+str(n_motif))

			if row[0][:6] == "#Motif":
				pos_neg = pos_neg * -1
			
			if row[0][:4] == ">chr":
				if pos_neg == 1:
					pos_motifs["Motif"+str(n_motif)] += 1
				elif pos_neg == -1:
					neg_motifs["Motif"+str(n_motif)] += 1			
		
		except IndexError:
			pass
	
	
	for k in keys:
		
		non_can_nkmers = non_can_seq_len_count  -  (kmer_len * non_can_seq_count)
		can_nkmers = can_seq_len_count  -  (kmer_len * can_seq_count)
		
		non_can_count = pos_motifs[k]
		can_count = neg_motifs[k]		
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
	main(sys.argv[1], sys.argv[2], sys.argv[3] ) 
