import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re
import random


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


def U2_score_calculator( E5, E3, U2_GTAG_3, U2_GTAG_5, U2_GTAG_5_max_score, U2_GTAG_3_max_score, TOTAL_U2_max_score):




	U2_score = 0
	ME5_U2_score = 0
	ME3_U2_score = 0	

	i = 0


	for N in E5:
		if N!="N":
			U2_score += U2_GTAG_3[N][i]
			ME5_U2_score += U2_GTAG_3[N][i]
			i += 1

	i = 0

	for N in E3:
		if N!="N":
			U2_score += U2_GTAG_5[N][i]
			ME3_U2_score += U2_GTAG_5[N][i]
			i += 1

	ME3_U2_score = percent(ME3_U2_score, U2_GTAG_5_max_score)
	ME5_U2_score = percent(ME5_U2_score, U2_GTAG_3_max_score)

	U2_score = percent(U2_score, TOTAL_U2_max_score)


	return ME3_U2_score, ME5_U2_score, U2_score


def main(bed12, U2_GTAG_5_file, U2_GTAG_3_file):

	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])
	
	TOTAL_U2_max_score = U2_GTAG_5_max_score + U2_GTAG_3_max_score


	for row in csv.reader(open(bed12), delimiter = '\t'):
		
		csv.field_size_limit(1000000000)

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))	


		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chrom = row[0]

		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = chrom + ":" +  str(istart) + row[5] + str(iend)

			intron_seq = Genome[chrom][istart:iend]

			if strand == '-':
				intron_seq = intron_seq.reverse_complement()

			intron_seq = str(intron_seq).upper()


			donor_matchs = []
			aceptor_matchs = []

			for donor in re.finditer("GT", intron_seq):

				donor_matchs.append(donor.start())


			for aceptor in re.finditer("AG", intron_seq):

				aceptor_matchs.append(aceptor.start())


			intron_SIM_MEs = []


			attempts = 10

			if len(aceptor_matchs) >= attempts: 

				for a in random.sample(aceptor_matchs, attempts):

					for ME_length in range(1,28):

						is_donor = a + 2 + ME_length

						if is_donor in set(donor_matchs):

							estart = a+2
							eend = is_donor

							E5 = intron_seq[estart-14:estart+3]
							E3 = intron_seq[eend-3:eend+10]


							SIM_ME = [ME_length, E5, E3]

							if len(E5)==17 and len(E3)==13:

								intron_SIM_MEs.append(SIM_ME)


			status="Simulated"


			if len(intron_SIM_MEs)>0:

				ME_length, E5, E3 = random.choice(intron_SIM_MEs)
				ME3_U2_score, ME5_U2_score, U2_score = U2_score_calculator(E5, E3, U2_GTAG_3, U2_GTAG_5, U2_GTAG_5_max_score, U2_GTAG_3_max_score, TOTAL_U2_max_score)

				print status, ME_length, E5, E3, ME5_U2_score, ME3_U2_score, U2_score


		###Geting the score for the annotated exons	

		unique_exons = set([])

		for q1, b in zip(qstarts[1:-1], blocksizes[1:-1]):
			estart = start + q1
			eend = start + q1 + b
			elenght = eend - estart

			E5 = Genome[chrom][estart-14:estart+3]
			E3 = Genome[chrom][eend-3:eend+10]

			if strand=="-":
				E5 = Genome[chrom][eend-3:eend+14].reverse_complement()
				E3 = Genome[chrom][estart-10:estart+3].reverse_complement()




			E5 = str(E5).upper()
			E3 = str(E3).upper()

			exon = "_".join(map(str, [chrom, strand, estart, eend]))


			status="Annotated"

			E3_U2_score, E5_U2_score, U2_score = U2_score_calculator(E5, E3, U2_GTAG_3, U2_GTAG_5, U2_GTAG_5_max_score, U2_GTAG_3_max_score, TOTAL_U2_max_score)

			exon_info = " ".join(map(str, [exon, status, elenght, E5, E3, E5_U2_score, E3_U2_score, U2_score]))

			unique_exons.add(exon_info)

		for e in unique_exons:

			exon, status, elenght, E5, E3, E5_U2_score, E3_U2_score, U2_score = e.split(" ")

			if set(E5)=={'A', 'C', 'G', 'T'} and set(E3)=={'A', 'C', 'G', 'T'}:

				print status, elenght, E5, E3, E5_U2_score, E3_U2_score, U2_score







if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4])


# python ~/my_src/PhD/ME/U2_simulation.py ../../../Genome/hg19/hg19.fa ../../../Genome/hg19/Tracks/Gene_annotation/gencode.v19.chr_patch_hapl_scaff.annotation.bed12 ../../../Genome/hg19/Tracks/SpliceRack/hg19_GT_AG_U2_5.good.matrix ../../../Genome/hg19/Tracks/SpliceRack/hg19_GT_AG_U2_3.good.matrix  > ME_U2.sim1.txt

#python ~/my_src/PhD/ME/U2_simulation.py $TEAM/Genome/mm10/mm10.fa $TEAM/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.annotation.bed12 $TEAM/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix $TEAM/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix  > ME_U2.sim1.txt



