import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import wWigIO
import numpy as np
from ngslib import IO
from scipy import stats
import re

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

def main(sim_fastq,  U2_GTAG_5_file, U2_GTAG_3_file, phylop_vertebrates, phylop_primates, exon_scores):

	MEs = set([])
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

	gencode_U2_scores = []
	gencode_mean_conservation_vertebrates = []
	gencode_mean_conservation_primates = []


	for row in csv.reader(open(exon_scores), delimiter = ' '):

		chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates = row

		gencode_U2_scores.append(float(U2_score))
		gencode_mean_conservation_vertebrates.append(float(mean_conservation_vertebrates))
		gencode_mean_conservation_primates.append(float(mean_conservation_primates))

	for row in csv.reader(open(sim_fastq), delimiter = '\t'):


		if row[0][0]=="@":

			SJ, ME_seq, estart, eend, total_coverage, n = row[0].split("_")

			len_ME = len(ME_seq)

			SJ = SJ[1:]
			SJ_chr, SJ_istart, SJ_iend = re.findall(r"[\w']+", SJ)


			SJ_len = int(SJ_iend) - int(SJ_istart)
			Kmer = SJ_len - (len_ME+1)
			P_ME = 1 - ( 1 - (float(1)/float(4**len_ME+4)))**Kmer	

			strand = "+"

			if "-" in SJ:
				strand = "-"

			estart = int(estart)
			eend = int(eend)

			MEs.add((SJ_chr, strand, estart, eend, P_ME))


	for m in MEs:

		chr, strand, estart, eend, P_ME = m

		estart, eend = sorted([estart, eend])

		E5 = str(Genome[chr][estart-14:estart+3]).upper()
		E3 = str(Genome[chr][eend-3:eend+10]).upper()


		if strand == "-":

			E5 = str(Genome[chr][eend-3:eend+14].reverse_complement()).upper()
			E3 = str(Genome[chr][estart-10:estart+3].reverse_complement()).upper()


		U2_score = 0

		i = 0

		for N in E5:
			if N!="N":
				U2_score += U2_GTAG_3[N][i]
				i += 1

		i = 0

		for N in E3:
			if N!="N":
				U2_score += U2_GTAG_5[N][i]
				i += 1

		U2_score = percent(U2_score, TOTAL_U2_max_score)

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

		ME_percentil_U2_score = stats.percentileofscore(gencode_U2_scores, U2_score)
		ME_percentil_mean_conservation_vertebrates = stats.percentileofscore(gencode_mean_conservation_vertebrates, mean_conservation_primates)
		ME_percentil_mean_conservation_primates = stats.percentileofscore(gencode_mean_conservation_primates, mean_conservation_vertebrates)

		overall_score = P_ME * (1- ME_percentil_U2_score/100) * (1 - ME_percentil_mean_conservation_vertebrates/100)

		if ME_percentil_mean_conservation_primates > ME_percentil_mean_conservation_vertebrates:
			overall_score = P_ME * (1- ME_percentil_U2_score/100) * (1 - ME_percentil_mean_conservation_primates/100)



		#print chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates
		print chr, estart, eend, strand, U2_score, ME_percentil_U2_score,  mean_conservation_vertebrates, ME_percentil_mean_conservation_vertebrates, mean_conservation_primates, ME_percentil_mean_conservation_primates, P_ME, overall_score

 #python ~/my_src/ME/Pipeline/percentile_scores_sim.py simulation_genome.fa sim_reads.fastq ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/hg19.100way.phyloP100way.bw ~/db/hg19.46way.phyloP46way.primates.bw

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main (sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]) 		