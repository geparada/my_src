import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from random import randint, sample
from operator import itemgetter

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0

def ascii_classifier(c): 

	ascii = ord(c)
	if ascii >= 48 and ascii <= 57:
		return 'number'
	else:
		return 'letter'

def cigar_parser(cigar):

	aux_str = ''
	cigar_vars = []

	for c in cigar:
		c_type = ascii_classifier(c)

		if c_type == 'number':
			aux_str += c

		elif c_type == 'letter':
			cigar_vars.append((c, int(aux_str)))
			aux_str = ''

	return cigar_vars

def main(sam):

	#micro_exons_found = defaultdict(set)
	total_micro_exons = defaultdict(set)

	micro_exons_found = set([])
	
	for row in csv.reader(open(sam), delimiter = '\t'):
		if row[0][0] != "@":
			
			read = row[0]
			flag = int(row[1])
			chr = row[2]
			start = int(row[3]) - 1             #Sam es 1 referenciado y es mas comodo trabajar en cordeneadas 0 refereciadas
			cigar = row[5]
			seq = row[9]
#			intron, micro_exon_seq, micro_exon_iend, micro_exon_istart, intron_coverage = read.split("_")
			intron, micro_exon_seq, micro_exon_iend, micro_exon_istart, intron_coverage, number_ID = read.split("_")  #en la nueva version de la simulacion se le agrego una variable




			micro_exon_length = len(micro_exon_seq)
			total_micro_exons[micro_exon_length].add(intron)

			alinments = []

			if (flag==0):

				alinments.append((read, flag, chr, start, cigar, seq))

#				if row[13].strip("AS:i:") == row[14].strip("XS:i:") and len(row)>=16 and "XA:Z:" in row[-1]:

#					c = 0

					#print row
					

#					for x in row[-1].strip("XA:Z:").split(";"):

#						if x!='':
#							tag, start, cigar, n = x.split(",")
#							if "+" in start:
#								start = int(start.strip("+")) - 1
								
#								alinments.append((read, flag, tag, start, cigar, seq))


			for i in alinments:

				read, flag, chr, start, cigar, seq = i
				intron, micro_exon_seq, micro_exon_iend, micro_exon_istart, intron_coverage, number_ID = read.split("_")


				cigar_list = cigar_parser(cigar)



				Exon_starts = [start]
				Exon_ends = []
				
				block = 0
				var_index = 0

				intron, micro_exon_seq, micro_exon_start, micro_exon_end, intron_coverage, n_intron_coverage = read.split("_")

				micro_exon_chr = intron.split(":")[0]


				micro_exon_start = int(micro_exon_start)
				micro_exon_end = int(micro_exon_end)

				micro_exon = (micro_exon_chr, micro_exon_start, micro_exon_end)


						
				for var in cigar_list:
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
						
					if var_index == len(cigar_list):
						Exon_ends.append(Exon_starts[-1] + block)

				### Extracting exonic sequences ###

				n = 0

				for block in cigar_list:

					block_type, block_len = block

					if block_type == "M":

						estart = start + n
						eend = start + n + block_len

						exon_found = (chr, estart, eend)

						if micro_exon == exon_found:

							micro_exons_found.add(micro_exon)


							#print read, block_type, estart, eend, micro_exon_chr, micro_exon_start, micro_exon_end

					n += block_len


	micro_exons_found_lengths = defaultdict(int)


	for i in micro_exons_found:

		micro_exon_chr, micro_exon_start, micro_exon_end = i

		micro_exons_found_lengths[micro_exon_end - micro_exon_start] += 1


	for l in range(25):

		print l+1, micro_exons_found_lengths[l+1], len(total_micro_exons[l+1]), percent(micro_exons_found_lengths[l+1], len(total_micro_exons[l+1]))
		



if __name__ == '__main__':
	main(sys.argv[1])