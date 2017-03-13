import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pyBigWig


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
	N_freq = []
	
	for row in reader:
		A = row[header_dict["A"]]
		C = row[header_dict["C"]]
		G = row[header_dict["G"]]
		T = row[header_dict["T"]]
		
		A_frec.append(float(A))
		C_frec.append(float(C))
		G_frec.append(float(G))
		T_frec.append(float(T))
		N_freq.append(0)
			
	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	matrix["N"] = N_freq
	
	return matrix


def main(ME_Final, U2_GTAG_5_file, U2_GTAG_3_file, phylop_vertebrates):

	reader = csv.reader(open(ME_Final), delimiter = '\t')

	headers = reader.next()

	phylop_vertebrates_bw = pyBigWig.open(phylop_vertebrates)

	phylop_vertebrates_bw = pyBigWig.open(phylop_vertebrates)


	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])


	for row in reader:


		ME, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_MEs, total_ME, ME_P_value, Filtered = row


		ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")

		ME_start = int(ME_start)
		ME_end = int(ME_end)


		L = 350

		up_seq = Genome[ME_chrom][ME_start-L:ME_start]
		down_seq = Genome[ME_chrom][ME_end:ME_end+L]

		if ME_strand=="-":

			
			up_seq = Genome[ME_chrom][ME_end:ME_end+L].reverse_complement()
			down_seq = Genome[ME_chrom][ME_start-L:ME_start].reverse_complement()



		up_seq = str(up_seq).upper()
		down_seq = str(down_seq).upper()


		# ME5 = str(Genome[ME_chr][ME_start-14:ME_start+3]).upper()
		# ME3 = str(Genome[ME_chr][ME_end-3:ME_end+10]).upper() 


		up_donor_U2_scores = []
		down_donor_U2_scores = []


		for i in range(L-13):

			up_donor = up_seq[i:i+13]
			down_donor = down_seq[i:i+13]


			up_donor_U2_score = 0
			down_donor_U2_score = 0

			p = 0

			for N in up_donor:
				up_donor_U2_score += U2_GTAG_5[N][p]
				p += 1

			p = 0

			for N in down_donor:
				down_donor_U2_score += U2_GTAG_5[N][p]
				p += 1


			up_donor_U2_score = percent( up_donor_U2_score, U2_GTAG_5_max_score)
			down_donor_U2_score = percent( down_donor_U2_score, U2_GTAG_5_max_score)

			up_donor_U2_scores.append(up_donor_U2_score)
			down_donor_U2_scores.append(down_donor_U2_score)


			up_donor_start = ME_start-L+i
			up_donor_end = ME_start-L+i+13

			down_donor_start = ME_end+i
			down_donor_end = ME_end+i+13

			up_donor_pylop = phylop_vertebrates_bw.stats(ME_chrom, up_donor_start, up_donor_end, type="mean")[0]
			down_donor_pylop = phylop_vertebrates_bw.stats(ME_chrom, down_donor_start, down_donor_end, type="mean")[0]

			if ME_strand=="-":

				up_donor_start = ME_end+L-i-13
				up_donor_end = ME_end+L-i

				down_donor_start = ME_start-i-13
				down_donor_end = ME_start-i

				up_donor_pylop = phylop_vertebrates_bw.stats(ME_chrom, up_donor_start, up_donor_end, type="mean")[0]
				down_donor_pylop = phylop_vertebrates_bw.stats(ME_chrom, down_donor_start, down_donor_end, type="mean")[0]



			if up_donor_pylop==None:
				up_donor_pylop=0

			if down_donor_pylop==None:
				down_donor_pylop=0


			if down_donor[3:5]=="GT":
				print ME, ME_chrom, ME_strand, i, down_donor_start, down_donor_end, down_donor,  down_donor_U2_score, down_donor_pylop

			#print ME_chrom, down_donor, Genome[ME_chrom][up_donor_pylop:down_donor_pylop], ME_st,  down_donor_U2_score, down_donor_pylop, 


		#print ME, up_donor_U2_scores, down_donor_U2_scores





if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

#python ~/my_src/PhD/ME/Stats/SS_scan.py /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/mm10.fa  ME_final_filtered.txt /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60way.bw 