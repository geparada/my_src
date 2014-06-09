import csv
import sys
import re
from collections import defaultdict
import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import scipy.stats


def fold (c, total):
	try:
		return (float(c))/float(total)
	except ZeroDivisionError:
		return 0

def chi_square (o, e):
	""" Calculate the chi-square """
	try:
		return ((o -e)**2)/e
	except ZeroDivisionError:
		return 0		


def main (non_can, can, fimo_non_can, fimo_can):
	""" Evalua significancia de los motivos encontrados """
	
	n_non_can = 0
	len_non_can = 0
	n_can = 0
	len_can = 0

	results = []

	reader1 = csv.reader(open(fimo_non_can), delimiter = '\t')
	reader2 = csv.reader(open(fimo_can), delimiter = '\t')	


	
	for record in SeqIO.parse(non_can, "fasta"):
		n_non_can += 1
		len_non_can += len(record.seq)
	
	for record in SeqIO.parse(can, "fasta"):
		n_can += 1
		len_can += len(record.seq)	
	
	reader1.next()
	reader2.next()
	
	
	motifs_len = {}
	
	motifs_count_non_can = defaultdict(int)	
	motifs_count_can = defaultdict(int)	
	
	
	for row in reader1:
		motif = row[0]
		seq = row[8]
		
		motifs_count_non_can[motif] += 1
		motifs_len[motif] = len(seq)

	for row in reader2:
		motif = row[0]
		seq = row[8]
		
		motifs_count_can[motif] += 1
		motifs_len[motif] = len(seq)
	
		
	for i in range(15):
		n = str(i+1)
		   
		if motifs_len.has_key(n)==False: #Para relenar el dicionario en caso de que el count del motif sea 0
			motifs_len[n] = 0
		
		k = motifs_len[n]
		nkmers_non_can = float(len_non_can - n_non_can * k) 
		nkmers_can = float(len_can - n_can * k)
		
		non_can_count = float(motifs_count_non_can[n])
		can_count = float(motifs_count_can[n])
		
		non_can_freq = non_can_count/nkmers_non_can
		can_freq = can_count/nkmers_can
		
		expected_non_can_count = (nkmers_non_can/nkmers_can) * can_count
		
		results.append((n, non_can_count, expected_non_can_count, chi_square(non_can_count, expected_non_can_count))) 
		
#		print n, non_can_count, expected_non_can_count, chi_square(non_can_count, expected_non_can_count), scipy.stats.chi2.sf(chi_square(non_can_count, expected_non_can_count), 14)

	GL = len(results) -1 
	
	for r in sorted(results, key=lambda a: a[3], reverse=True):
		cluster = r[0]
		non_can = r[1]
		expected_non_can = r[2]
		chi2 = r[3]
		p_value = scipy.stats.chi2.sf(chi2, GL)
		
		print cluster, non_can, expected_non_can, chi2, p_value


		


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])	
