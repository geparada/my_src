import sys
import csv
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
import re
import wWigIO
from ngslib import BigWigFile
import numpy as np
from scipy import stats


Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0

def bigwig_fetcher(bigwig, ichr, istart, iend):

	bw=BigWigFile(bigwig)
	
	scores = []

	with BigWigFile(bigwig) as bw:

		for i in bw.fetch(chrom=ichr,start=istart,stop=iend):

			scores.append(i.score)

		return np.mean(scores)

	bw.wWigIO.close()				

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


def SJ_score(max_score, PWM_5, PWM_3, SEQ5, SEQ3):

	score = 0

	i = 0

	for N in SEQ5:
		score += PWM_3[N][i]
		i += 1

	i = 0

	for N in SEQ3:
		score += PWM_5[N][i]
		i += 1

	score = percent(score, max_score)

	return score

def main(gencode_bed, U2_GTAG_5_file, U2_GTAG_3_file, U2_GCAG_5_file, U2_GCAG_3_file, U12_ATAC_5_file, U12_ATAC_3_file, phylop_vertebrates, phylop_primates, simulation_scores):

	donors = defaultdict(set)
	acepors = defaultdict(set)

	putative_ME = set([])

	simulation_U2 = []
	simulation_vertebrates = []
	simulation_primates = []

	for row in csv.reader(open(simulation_scores), delimiter = ' '):

		chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates = row

		simulation_U2.append(float(U2_score))
		simulation_vertebrates.append(float(mean_conservation_vertebrates))
		simulation_primates.append(float(mean_conservation_primates))

	for row in csv.reader(open(gencode_bed), delimiter = '\t'):
		

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chrom = row[0]
		qstart = 0

		for q1, q2, b, b2 in zip(qstarts, qstarts[1:], blocksizes, blocksizes[1:]):
			
			qstart = qstart + b
			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = chrom + ":" + str(istart) + strand + str(iend)
			ilength = iend - istart

			info = " ".join(map(str, [chrom, istart, iend, strand]))

			donor_ID = " ".join(map(str, [chrom, istart, strand]))
			aceptor_ID = " ".join(map(str, [chrom, iend, strand]))

									 
			if ilength >= 80:  # hay que agregarle el filtro de los micro exones!!

				donors[donor_ID].add(iend)
				acepors[aceptor_ID].add(istart)


	ME_matches = defaultdict(int)
	alt_introns = {}


	for d in donors.items():

		donor, iends = d

		chrom, istart, strand = donor.split(" ")
		istart = int(istart)

		for ME_start in iends:

			for ME_end in iends:

				if 0 < ME_end - ME_start <= 25:

					iend = ME_start

					intron = chrom + ":" + str(istart) + strand + str(iend)
					alt_intron = chrom + ":" + str(istart) + strand + str(ME_end)

					ME_seq = str(Genome[chrom][ME_start:ME_end]).upper()
					intron_seq = str(Genome[chrom][istart:iend]).upper()

					if strand == "-":
						ME_seq = str(Genome[chrom][ME_start:ME_end].reverse_complement()).upper()
						intron_seq = str(Genome[chrom][istart:iend].reverse_complement()).upper()

					seq = str("AG" + ME_seq + "GT").upper()

					ME_info = ""


					c = 0

					while intron_seq.find(seq, c)!=-1:

						c = intron_seq.find(seq, c)
						n = c + 2
							
						ME_info = " ".join([intron, ME_seq, chrom, str(istart+n), strand, str(istart+n+len(ME_seq)), "alternative_aceptor"])
						
						if strand == "-":
							 ME_info = " ".join([intron, ME_seq, chrom, str(iend-n-len(ME_seq)), strand, str(iend-n), "alternative_donor"])
						

						putative_ME.add(ME_info)
						alt_introns[intron] = alt_intron
						ME_matches[intron] += 1


						c += 1



	for a in acepors.items():

		aceptor, istarts = a

		chrom, iend, strand = aceptor.split(" ")
		iend = int(iend)



		for ME_start in istarts:

			for ME_end in istarts:

				if 0 < ME_end - ME_start <= 25:

					istart = ME_end

					intron = chrom + ":" + str(istart) + strand + str(iend)
					alt_intron = chrom + ":" + str(ME_start) + strand + str(iend)

					ME_seq = str(Genome[chrom][ME_start:ME_end]).upper()
					intron_seq = str(Genome[chrom][istart:iend]).upper()

					if strand == "-":
						ME_seq = str(Genome[chrom][ME_start:ME_end].reverse_complement()).upper()
						intron_seq = str(Genome[chrom][istart:iend].reverse_complement()).upper()

					seq = str("AG" + ME_seq + "GT").upper()

					ME_info = ""

					c = 0

					while intron_seq.find(seq, c)!=-1:

						c = intron_seq.find(seq, c)
						n = c + 2
							
						ME_info = " ".join([intron, ME_seq, chrom, str(istart+n), strand, str(istart+n+len(ME_seq)), "alternative_donor"])
						
						if strand == "-":
							 ME_info = " ".join([intron, ME_seq, chrom, str(iend-n-len(ME_seq)), strand, str(iend-n), "alternative_aceptor"])
						

						putative_ME.add(ME_info)
						alt_introns[intron] = alt_intron
						ME_matches[intron] += 1


						c += 1

	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])
	
	TOTAL_GTAG_max_score = U2_GTAG_5_max_score + U2_GTAG_3_max_score

	U2_GCAG_5 = PWM_to_dict(U2_GCAG_5_file)
	U2_GCAG_3 = PWM_to_dict(U2_GCAG_3_file)

	U2_GCAG_5_max_score = 0
	U2_GCAG_3_max_score = 0

	for index in range(13):
		U2_GCAG_5_max_score += max(U2_GCAG_5['A'][index], U2_GCAG_5['C'][index], U2_GCAG_5['T'][index], U2_GCAG_5['G'][index])

	for index in range(17):
		U2_GCAG_3_max_score += max(U2_GCAG_3['A'][index], U2_GCAG_3['C'][index], U2_GCAG_3['T'][index], U2_GCAG_3['G'][index])
	
	TOTAL_GCAG_max_score = U2_GCAG_5_max_score + U2_GCAG_3_max_score

	U12_ATAC_5 = PWM_to_dict(U12_ATAC_5_file)
	U12_ATAC_3 = PWM_to_dict(U12_ATAC_3_file)

	U12_ATAC_5_max_score = 0
	U12_ATAC_3_max_score = 0

	for index in range(13):
		U12_ATAC_5_max_score += max(U12_ATAC_5['A'][index], U12_ATAC_5['C'][index], U12_ATAC_5['T'][index], U12_ATAC_5['G'][index])

	for index in range(17):
		U12_ATAC_3_max_score += max(U12_ATAC_3['A'][index], U12_ATAC_3['C'][index], U12_ATAC_3['T'][index], U12_ATAC_3['G'][index])
	
	TOTAL_ATAC_max_score = U12_ATAC_5_max_score + U12_ATAC_3_max_score

	for e in putative_ME:


		intron, ME_seq, chrom, ME_start, strand, ME_end, ME_type = e.split(" ")


		chrom, istart, iend = re.findall(r"[\w']+", intron)

		istart, iend, ME_start, ME_end = map(int, [istart, iend, ME_start, ME_end])

		alt_intron = alt_introns[intron]
		total_number_of_micro_exons_matches = ME_matches[intron]

		dn = str(Genome[chrom][istart:(istart+2)] + Genome[chrom][(iend-2):iend]).upper()
		ME5 = str(Genome[chrom][ME_start-14:ME_start+3]).upper()
		ME3 = str(Genome[chrom][ME_end-3:ME_end+10]).upper()

		I3 = str(Genome[chrom][iend-14:iend+3]).upper()
		I5 = str(Genome[chrom][istart-3:istart+10]).upper()


		mean_conservation_vertebrates = bigwig_fetcher(phylop_vertebrates, chrom, ME_start, ME_end)
		mean_conservation_primates = bigwig_fetcher(phylop_primates, chrom, ME_start, ME_end)


		if strand == "-":

			ME5 = str(Genome[chrom][ME_end-3:ME_end+14].reverse_complement()).upper()
			ME3 = str(Genome[chrom][ME_start-10:ME_start+3].reverse_complement()).upper()
			dn = str((Genome[chrom][istart:(istart+2)] + Genome[chrom][(iend-2):iend]).reverse_complement()).upper()

			I3 = str(Genome[chrom][istart-3:istart+14].reverse_complement()).upper()
			I5 = str(Genome[chrom][iend-10:iend+3].reverse_complement()).upper()

		ME_SJ_score = SJ_score(TOTAL_GTAG_max_score, U2_GTAG_5, U2_GTAG_3, ME5, ME3)

		I_GTAG_SJ_score = SJ_score(TOTAL_GTAG_max_score, U2_GTAG_5, U2_GTAG_3, I3, I5)
		I_GCAG_SJ_score = SJ_score(TOTAL_GCAG_max_score, U2_GCAG_5, U2_GCAG_3, I3, I5)
		I_ATAC_SJ_score = SJ_score(TOTAL_ATAC_max_score, U12_ATAC_5, U2_GTAG_3, I3, I5)

		I_SJ_score = max([I_GTAG_SJ_score, I_GCAG_SJ_score, I_ATAC_SJ_score])

		#if dn!="GTAG" and dn!="GCAG" and dn!="ATAC" and ME_SJ_score>=80:

		U2_FPR = 1 - stats.percentileofscore(simulation_U2, float(ME_SJ_score))/100
		vertebrates_FPR = 1 - stats.percentileofscore(simulation_vertebrates, float(mean_conservation_vertebrates))/100
		primates_FPR = 1 - stats.percentileofscore(simulation_primates, float(mean_conservation_primates))/100

		total_score = U2_FPR * vertebrates_FPR * float(total_number_of_micro_exons_matches)

		if 	ME_SJ_score > I_SJ_score:

			#print e, dn, mean_conservation_vertebrates, mean_conservation_primates, ME_SJ_score, I_SJ_score, total_number_of_micro_exons_matches, total_score 

			if total_score <= 0.05 and U2_FPR <= 0.2:
				#print e, dn,mean_conservation_vertebrates, mean_conservation_primates, ME_SJ_score, I_SJ_score, alt_intron, total_number_of_micro_exons_matches, total_score 
				if len(ME_seq) >= 7:
					print e, dn,mean_conservation_vertebrates, mean_conservation_primates, ME_SJ_score, I_SJ_score, alt_intron, total_number_of_micro_exons_matches, total_score 
				elif total_score <= 0.005 and vertebrates_FPR <= 0.05 and U2_FPR <= 0.2:
					print e, dn,mean_conservation_vertebrates, mean_conservation_primates, ME_SJ_score, I_SJ_score, alt_intron, total_number_of_micro_exons_matches, total_score 

			
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])	


#python  ~/my_src/ME/Pipeline/Round_1/miss_annotated_exon_gencode.py ~/db/genome/hg19.fa ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/PWM/hg19_GC_AG_U2_5.good.matrix ~/db/PWM/hg19_GC_AG_U2_3.good.matrix ~/db/PWM/hg19_AT_AC_U12_5.good.matrix ~/db/PWM/hg19_AT_AC_U12_3.good.matrix ~/db/hg19.100way.phyloP100way.bw ~/db/hg19.46way.phyloP46way.primates.bw /media/HD3/Resultados/Micro_exons/simulation/ROC/ME_sim.scores 
