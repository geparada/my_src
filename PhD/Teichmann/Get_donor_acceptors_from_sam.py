import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict


csv.field_size_limit(500 * 1024 * 1024)


def ascii_classifier(c): 

	ascii = ord(c)
	if ascii >= 48 and ascii <= 57:
		return 'number'
	else:
		return 'letter'
			

def main(sam, anchor):          #hay que indicar s es Rd1 o Rd2
	reader = csv.reader(open(sam), delimiter = '\t')
	
	
	SS = defaultdict(set)

	for row in reader:
		if row[0][0] != "@":
			if "N" in row[5]:
				read = row[0]
				flag = row[1]
				chr = row[2]
				start = int(row[3]) - 1             #Sam es 1 referenciado y es mas comodo trabajar en cordeneadas 0 refereciadas
				cigar = row[5]
				seq = row[9]
						
				aux_str = ''
				cigar_vars = [] 
							
				for c in cigar:
					c_type = ascii_classifier(c)

					if c_type == 'number':
						aux_str += c

					elif c_type == 'letter':
						cigar_vars.append((c, int(aux_str)))
						aux_str = ''

				Exon_starts = [start]
				Exon_ends = []
				
				block = 0
				var_index = 0
						
				for var in cigar_vars:
					var_type = var[0]
					var_value = var[1]
					var_index += 1
					
					if var_type == 'M':
						block += var_value						
						
					if var_type == 'D':
						block += var_value
											
					if var_type == 'I':
						block += 0
						
					if var_type == 'N':
						Exon_ends.append(Exon_starts[-1] + block)
						Exon_starts.append(Exon_ends[-1] + var_value)
						block = 0
						
					if var_index == len(cigar_vars):
						Exon_ends.append(Exon_starts[-1] + block)

				
				for e5s, e5e, e3s, e3e  in zip(Exon_starts, Exon_ends, Exon_starts[1:], Exon_ends[1:]):
					e5len= e5e - e5s 
					e3len = e3e - e3s
					istart = e5e
					iend = e3s
					ilen = iend - istart
					intron = "_".join([chr, str(istart), str(iend)])


					if e5len >= anchor <= e3len:

						SS["_".join([chr, str(istart)])].add(seq)
						SS["_".join([chr, str(iend)])].add(seq)
						

	for i in SS.items():
		ss, seqs = i
		cov = len(seqs)

		print ss, cov

						
			



if __name__ == '__main__':
	main(sys.argv[1],  int(sys.argv[2]))
