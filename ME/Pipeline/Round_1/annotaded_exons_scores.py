import sys
import csv
from collections import defaultdict
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import wWigIO
import numpy as np
from ngslib import IO



Genome = {}

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0		


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

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

def PWM_to_dict(file):
	reader = csv.reader(open(file), delimiter = '\t')
	header = reader.next()
	header_dict = {}
	col = 0
	
	matrix = {}
	
	for name in header:
		header_dict[name] = col
		col += 1
	
	A_frec = []
	C_frec = []
	G_frec = []
	T_frec = []
	
	for row in reader:
		A = row[header_dict["A"]]
		C = row[header_dict["C"]]
		G = row[header_dict["G"]]
		T = row[header_dict["T"]]
		
		A_frec.append(float(A))
		C_frec.append(float(C))
		G_frec.append(float(G))
		T_frec.append(float(T))
			
	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	
	return matrix


def main(gencode_bed, U2_GTAG_5_file, U2_GTAG_3_file, phylop_vertebrates, phylop_primates):

	wWigIO.open(phylop_vertebrates)
	wWigIO.open(phylop_primates)

	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])
	
	TOTAL_U2_max_score = U2_GTAG_5_max_score + U2_GTAG_3_max_score

	exons = set([])

	for row in csv.reader(open(gencode_bed), delimiter = '\t'):
		
		csv.field_size_limit(1000000000)

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]

		

		for q1, b in zip(qstarts[1:-1], blocksizes[1:-1]):
			estart = start + q1
			eend = start + q1 + b


			E5 = str(Genome[chr][estart-14:estart+3]).upper()
			E3 = str(Genome[chr][eend-3:eend+10]).upper()


			if strand == "-":

				E5 = str(Genome[chr][eend-3:eend+14].reverse_complement()).upper()
				E3 = str(Genome[chr][estart-10:estart+3].reverse_complement()).upper()


			U2_score = 0

			i = 0

			for N in E5:
				U2_score += U2_GTAG_3[N][i]
				i += 1

			i = 0

			for N in E3:
				U2_score += U2_GTAG_5[N][i]
				i += 1

			U2_score = percent(U2_score, TOTAL_U2_max_score)



			if E5[-5:-3]=="AG" and E3[3:5] == "GT":



				exons.add((chr, estart, eend, strand, U2_score))

			# if " ".join([chr, estart, eend]) == "chr17 26597935 26598725":
			# 	print 


	for e in exons:



		chr, estart, eend, strand, U2_score = e

		conservation_vertebrates = wWigIO.getIntervals(phylop_vertebrates, chr, estart-2, eend+2)
		conservation_primates = wWigIO.getIntervals(phylop_primates, chr, estart-2, eend+2)

		mean_conservation_vertebrates = 0
		mean_conservation_primates = 0

		for i in conservation_vertebrates:

			mean_conservation_vertebrates += i[2]

		try:

			mean_conservation_vertebrates = mean_conservation_vertebrates/len(conservation_vertebrates)

		except ZeroDivisionError:
			pass

		
		for i in conservation_primates:

			mean_conservation_primates += i[2]

		try:

			mean_conservation_primates = mean_conservation_primates/len(conservation_primates)

		except ZeroDivisionError:
			pass

		print chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])


#python ~/my_src/ME/Pipeline/annotaded_exons_scores.py ~/db/genome/hg19.fa ~/db/transcriptome/hg19/Gene_models/gencode/v17/gencode.v17.annotation.bed12 ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/hg19.100way.phyloP100way.bw ~/db/hg19.46way.phyloP46way.primates.bw