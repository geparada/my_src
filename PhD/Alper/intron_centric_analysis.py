import csv
import sys
from collections import defaultdict
import numpy as np
import HTSeq
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re

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

def PWM_score(seq, PWM_dict):

	PWM_max_score = 0
	index = 0
	for N in range(len(PWM_dict["A"])):
		PWM_max_score += max(PWM_dict['A'][index], PWM_dict['C'][index], PWM_dict['T'][index], PWM_dict['G'][index])
		index += 1


	PWM_score = 0
	index = 0
	for N in seq:
		try:
			frec = PWM_dict[N][index]
			PWM_score += frec
			index += 1
		except KeyError:
			pass

	percentual_PWM_score = percent(PWM_score, PWM_max_score)

	return percentual_PWM_score




def rpkm_calculator(counts, length, lib_size):

	counts += 1 #pseudo count

	rpkm = (counts * 1000)/(length*lib_size)

	return rpkm


def z_score(x, mean, std):

	return float(x - mean)/float(std + 0.00000001)


def main(intronic_counts, exonic_counts, EISA, ce_RepeatMasker, intronic_gtf, exonic_gtf, DEXSeq_out, U2_GTAG_5_file, U2_GTAG_3_file):


	#Only the introns and exons present on the GTF files have to be considered

	gene_introns = defaultdict(list)
	gene_exons =  defaultdict(list)


	for row in csv.reader(open(intronic_gtf), delimiter = '\t'):


		chrom, annotation, feature, start, end, score, strand, cero, ID = row
		gene = ID.split(" ")[1].strip(";")
		intron = "_".join([chrom, start, end, strand])

		if feature == "exon":
			gene_introns[gene].append(intron)

	for row in csv.reader(open(exonic_gtf), delimiter = '\t'):

		chrom, annotation, feature, start, end, score, strand, cero, ID = row
		gene = ID.split(" ")[1].strip(";")
		exon = "_".join([chrom, start, end, strand])

		if feature == "exon":
			gene_exons[gene].append(exon)


	reader_intronic1 = csv.reader(open(intronic_counts), delimiter = '\t')
	reader_exonic1 = csv.reader(open(exonic_counts), delimiter = '\t')
	reader_intronic2 = csv.reader(open(intronic_counts), delimiter = '\t')
	reader_exonic2 = csv.reader(open(exonic_counts), delimiter = '\t')			
	reader_EISA = csv.reader(open(EISA), delimiter = ' ')
	reader_RepeatMasker = csv.reader(open(ce_RepeatMasker), delimiter = '\t')
	reader_DEXSeq = csv.reader(open(DEXSeq_out), delimiter = ' ')

	reader_intronic1.next()
	reader_exonic1.next()
	reader_intronic2.next()      
	reader_exonic2.next()

	reader_intronic1.next()
	reader_exonic1.next()
	reader_intronic2.next()      
	reader_exonic2.next()

	reader_EISA.next()
	reader_RepeatMasker.next()
	reader_DEXSeq.next()

	exonic_EISA = []
	intronic_EISA = []

	genes_exonic_EISA = set([])
	genes_intronic_EISA = set([])

	for row in reader_EISA:

		gene, Ex_WT_1, Ex_WT_2, Ex_WT_3, Ex_emb4a_1, Ex_emb4a_2, Ex_emb4a_3, Ex_emb4b_1, Ex_emb4b_2, Ex_emb4b_3, In_WT_1, In_WT_2, In_WT_3, In_emb4a_1, In_emb4a_2, In_emb4a_3, In_emb4b_1, In_emb4b_2, In_emb4b_3, D_ex, D_in, EISA_logFC, logCPM, LR, PValue, EISA_FDR, baseMean, DESeq_log2FoldChange, lfcSE, stat, pvalue, DESeq_padj =  row

		EISA_FDR = float(EISA_FDR)
		DESeq_padj = float(DESeq_padj)

		EISA_logFC = float(EISA_logFC)
		DESeq_log2FoldChange = float(DESeq_log2FoldChange)



		if EISA_FDR < 0.05:

			if EISA_logFC > 0:

				exonic_EISA.append([gene, EISA_logFC, EISA_FDR, DESeq_log2FoldChange, EISA_logFC])
				genes_exonic_EISA.add(gene)

			else:

				intronic_EISA.append([gene, EISA_logFC, EISA_FDR, DESeq_log2FoldChange, EISA_logFC])
				genes_intronic_EISA.add(gene)


	lib_sizes = defaultdict(int)

	for row in reader_intronic1:

		ID, chrom, start, end, strand, length, SX1316_1, SX1316_2, SX1316_3, SX2000_1, SX2000_2, SX2000_3, SX2929_1, SX2929_2, SX2929_3, SX2930_1, SX2930_2, SX2930_3 = row

		SX1316_1 = float(SX1316_1)		
		SX1316_2 = float(SX1316_2)
		SX1316_3 = float(SX1316_3)

		SX2000_1 = float(SX2000_1)		
		SX2000_2 = float(SX2000_2)
		SX2000_3 = float(SX2000_3)

		SX2929_1 = float(SX2929_1)		
		SX2929_2 = float(SX2929_2)
		SX2929_3 = float(SX2929_3)

		SX2930_1 = float(SX2930_1)		
		SX2930_2 = float(SX2930_2)
		SX2930_3 = float(SX2930_3)		

		lib_sizes["SX1316_1"] += SX1316_1		
		lib_sizes["SX1316_2"] += SX1316_2
		lib_sizes["SX1316_3"] += SX1316_3

		lib_sizes["SX2000_1"] += SX2000_1		
		lib_sizes["SX2000_2"] += SX2000_2
		lib_sizes["SX2000_3"] += SX2000_3

		lib_sizes["SX2929_1"] += SX2929_1		
		lib_sizes["SX2929_2"] += SX2929_2
		lib_sizes["SX2929_3"] += SX2929_3

		lib_sizes["SX2930_1"] += SX2930_1		
		lib_sizes["SX2930_2"] += SX2930_2
		lib_sizes["SX2930_3"] += SX2930_3


	for row in reader_exonic1:

		ID, chrom, start, end, strand, length, SX1316_1, SX1316_2, SX1316_3, SX2000_1, SX2000_2, SX2000_3, SX2929_1, SX2929_2, SX2929_3, SX2930_1, SX2930_2, SX2930_3 = row

		SX1316_1 = float(SX1316_1)		
		SX1316_2 = float(SX1316_2)
		SX1316_3 = float(SX1316_3)

		SX2000_1 = float(SX2000_1)		
		SX2000_2 = float(SX2000_2)
		SX2000_3 = float(SX2000_3)

		SX2929_1 = float(SX2929_1)		
		SX2929_2 = float(SX2929_2)
		SX2929_3 = float(SX2929_3)

		SX2930_1 = float(SX2930_1)		
		SX2930_2 = float(SX2930_2)
		SX2930_3 = float(SX2930_3)		

		lib_sizes["SX1316_1"] += SX1316_1		
		lib_sizes["SX1316_2"] += SX1316_2
		lib_sizes["SX1316_3"] += SX1316_3

		lib_sizes["SX2000_1"] += SX2000_1		
		lib_sizes["SX2000_2"] += SX2000_2
		lib_sizes["SX2000_3"] += SX2000_3

		lib_sizes["SX2929_1"] += SX2929_1		
		lib_sizes["SX2929_2"] += SX2929_2
		lib_sizes["SX2929_3"] += SX2929_3

		lib_sizes["SX2930_1"] += SX2930_1		
		lib_sizes["SX2930_2"] += SX2930_2
		lib_sizes["SX2930_3"] += SX2930_3

	intronic_log2_emb4 = {}
	intron_ID = {}	

	for row in reader_intronic2:

		ID, chrom, start, end, strand, length, SX1316_1, SX1316_2, SX1316_3, SX2000_1, SX2000_2, SX2000_3, SX2929_1, SX2929_2, SX2929_3, SX2930_1, SX2930_2, SX2930_3 = row

		SX1316_1 = float(SX1316_1)		
		SX1316_2 = float(SX1316_2)
		SX1316_3 = float(SX1316_3)

		SX2000_1 = float(SX2000_1)		
		SX2000_2 = float(SX2000_2)
		SX2000_3 = float(SX2000_3)

		SX2929_1 = float(SX2929_1)		
		SX2929_2 = float(SX2929_2)
		SX2929_3 = float(SX2929_3)

		SX2930_1 = float(SX2930_1)		
		SX2930_2 = float(SX2930_2)
		SX2930_3 = float(SX2930_3)

		length = float(length)


		SX1316_1_rpkm = rpkm_calculator(SX1316_1, length, lib_sizes["SX1316_1"])
		SX1316_2_rpkm = rpkm_calculator(SX1316_2, length, lib_sizes["SX1316_2"])		
		SX1316_3_rpkm = rpkm_calculator(SX1316_3, length, lib_sizes["SX1316_3"])

		SX2000_1_rpkm = rpkm_calculator(SX2000_1, length, lib_sizes["SX2000_1"])
		SX2000_2_rpkm = rpkm_calculator(SX2000_2, length, lib_sizes["SX2000_2"])		
		SX2000_3_rpkm = rpkm_calculator(SX2000_3, length, lib_sizes["SX2000_3"])

		SX2929_1_rpkm = rpkm_calculator(SX2929_1, length, lib_sizes["SX2929_1"])
		SX2929_2_rpkm = rpkm_calculator(SX2929_2, length, lib_sizes["SX2929_2"])		
		SX2929_3_rpkm = rpkm_calculator(SX2929_3, length, lib_sizes["SX2929_3"])

		SX2930_1_rpkm = rpkm_calculator(SX2930_1, length, lib_sizes["SX2930_1"])
		SX2930_2_rpkm = rpkm_calculator(SX2930_2, length, lib_sizes["SX2930_2"])		
		SX2930_3_rpkm = rpkm_calculator(SX2930_3, length, lib_sizes["SX2930_3"])

		SX1316_rpkm_mean =  np.mean([SX1316_1_rpkm, SX1316_2_rpkm, SX1316_3_rpkm])
		SX2000_rpkm_mean =  np.mean([SX2000_1_rpkm, SX2000_2_rpkm, SX2000_3_rpkm])
		SX2929_rpkm_mean =  np.mean([SX2929_1_rpkm, SX2929_2_rpkm, SX2929_3_rpkm])
		SX2930_rpkm_mean =  np.mean([SX2930_1_rpkm, SX2930_2_rpkm, SX2930_3_rpkm])
		SX2929_SX2930_rpkm_mean =  np.mean([SX2929_1_rpkm, SX2929_2_rpkm, SX2929_3_rpkm, SX2930_1_rpkm, SX2930_2_rpkm, SX2930_3_rpkm])

		log2_emb4 =  np.log2(SX2929_SX2930_rpkm_mean / SX1316_rpkm_mean )

		intron = "_".join([chrom, start, end, strand])

		intronic_log2_emb4[intron] = log2_emb4

		intron_ID[ID] = intron		

	exonic_log2_emb4 = {}	

	for row in reader_exonic2:

		ID, chrom, start, end, strand, length, SX1316_1, SX1316_2, SX1316_3, SX2000_1, SX2000_2, SX2000_3, SX2929_1, SX2929_2, SX2929_3, SX2930_1, SX2930_2, SX2930_3 = row

		SX1316_1 = float(SX1316_1)		
		SX1316_2 = float(SX1316_2)
		SX1316_3 = float(SX1316_3)

		SX2000_1 = float(SX2000_1)		
		SX2000_2 = float(SX2000_2)
		SX2000_3 = float(SX2000_3)

		SX2929_1 = float(SX2929_1)		
		SX2929_2 = float(SX2929_2)
		SX2929_3 = float(SX2929_3)

		SX2930_1 = float(SX2930_1)		
		SX2930_2 = float(SX2930_2)
		SX2930_3 = float(SX2930_3)

		length = float(length)


		SX1316_1_rpkm = rpkm_calculator(SX1316_1, length, lib_sizes["SX1316_1"])
		SX1316_2_rpkm = rpkm_calculator(SX1316_2, length, lib_sizes["SX1316_2"])		
		SX1316_3_rpkm = rpkm_calculator(SX1316_3, length, lib_sizes["SX1316_3"])

		SX2000_1_rpkm = rpkm_calculator(SX2000_1, length, lib_sizes["SX2000_1"])
		SX2000_2_rpkm = rpkm_calculator(SX2000_2, length, lib_sizes["SX2000_2"])		
		SX2000_3_rpkm = rpkm_calculator(SX2000_3, length, lib_sizes["SX2000_3"])

		SX2929_1_rpkm = rpkm_calculator(SX2929_1, length, lib_sizes["SX2929_1"])
		SX2929_2_rpkm = rpkm_calculator(SX2929_2, length, lib_sizes["SX2929_2"])		
		SX2929_3_rpkm = rpkm_calculator(SX2929_3, length, lib_sizes["SX2929_3"])

		SX2930_1_rpkm = rpkm_calculator(SX2930_1, length, lib_sizes["SX2930_1"])
		SX2930_2_rpkm = rpkm_calculator(SX2930_2, length, lib_sizes["SX2930_2"])		
		SX2930_3_rpkm = rpkm_calculator(SX2930_3, length, lib_sizes["SX2930_3"])

		SX1316_rpkm_mean =  np.mean([SX1316_1_rpkm, SX1316_2_rpkm, SX1316_3_rpkm])
		SX2000_rpkm_mean =  np.mean([SX2000_1_rpkm, SX2000_2_rpkm, SX2000_3_rpkm])
		SX2929_rpkm_mean =  np.mean([SX2929_1_rpkm, SX2929_2_rpkm, SX2929_3_rpkm])
		SX2930_rpkm_mean =  np.mean([SX2930_1_rpkm, SX2930_2_rpkm, SX2930_3_rpkm])
		SX2929_SX2930_rpkm_mean =  np.mean([SX2929_1_rpkm, SX2929_2_rpkm, SX2929_3_rpkm, SX2930_1_rpkm, SX2930_2_rpkm, SX2930_3_rpkm])

		log2_emb4 =  np.log2(SX2929_SX2930_rpkm_mean / SX1316_rpkm_mean )

		exon = "_".join([chrom, start, end, strand])

		exonic_log2_emb4[exon] = log2_emb4	

##### Are there any overepresented repeat over the downregulated introns belonging to the exonic EISA dots??


	# min_intron_log2_emb4_info = []

	# for i in gene_introns.items():

	# 	gene, introns = i

	# 	if gene in genes_exonic_EISA:

	# 		introns_log2_emb4 = []

	# 		for intron in introns:

	# 			log2_emb4 = intronic_log2_emb4[intron]
	# 			introns_log2_emb4.append((intron, log2_emb4))



	# 		min_intron_log2_emb4 = min(introns_log2_emb4, key = lambda x: x[1])[0]




	# 		chrom, start, end, strand = min_intron_log2_emb4.split("_")

	# 		min_intron_log2_emb4_interval = HTSeq.GenomicInterval( chrom, int(start), int(end), "." )

	# 		min_intron_log2_emb4_info.append([chrom, start, end, strand, min_intron_log2_emb4_interval, gene])



	# rep_min_intron = set([])
	# rep_min_intron_count = defaultdict(int)
						




	# for row in reader_RepeatMasker:


	# 	bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, rep_strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id = row
	# 	rep_interval = HTSeq.GenomicInterval( genoName, int(genoStart), int(genoEnd), "." )


	# 	for i in min_intron_log2_emb4_info:

	# 		chrom, start, end, strand, min_intron_log2_emb4_interval, gene = i

	# 		if min_intron_log2_emb4_interval.overlaps(rep_interval):

	# 			rep_min_intron.add(" ".join([gene, min_intron_log2_emb4, repName, repClass, repFamily]))




	# for i in rep_min_intron:


	# 	gene, min_intron_log2_emb4, repName, repClass, repFamily = i.split(" ")

	# 	rep_min_intron_count[repName] += 1


	# rep_min_intron_count = rep_min_intron_count.items()
	# rep_min_intron_count.sort(key=lambda x: x[1])

	# for i in reversed(rep_min_intron_count):

	# 	print "\t".join(map(str, i))


	#### Intron retention analysis


	DEXSeq_introns = set([])


	for row in reader_DEXSeq:



		groupID, featureID, exonBaseMean, dispersion, stat, pvalue, padj, control, knockout, log2fold_control_knockout, genomicData_seqnames, genomicData_start, genomicData_end, genomicData_width, genomicData_strand, countData_mutant1, countData_untreated1, countData_untreated2, countData_mutant2, countData_untreated3, countData_mutant3, transcripts = row


		ID =  groupID + featureID.replace("E", "_")

		DEXSeq_introns.add(intron_ID[ID])


	intron_CELE45 = {}


	for row in reader_RepeatMasker:


		bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, rep_strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id = row
		
		if repName=="CELE45":

			rep_interval = HTSeq.GenomicInterval( genoName, int(genoStart), int(genoEnd), "." )

			for g in gene_introns.items():

					gene, introns = g

					for intron in introns:

						try:

							chrom, start, end, strand =  intron.split("_")
							intron_interval = HTSeq.GenomicInterval( chrom, int(start), int(end), "." )

						except ValueError:
							pass

						if intron_interval.overlaps(rep_interval):

							intron_CELE45[intron] = milliDiv




	gene_mean_std_exon_log2_emb4 = {}
	

	for i in gene_exons.items():

		gene, exons = i

		exons_log2_emb4 = []

		n = 0

		for exon in exons:

			n+=1

			log2_emb4 = exonic_log2_emb4[exon]
			exons_log2_emb4.append(log2_emb4)

		mean_exons_log2_emb4 = np.mean(exons_log2_emb4)
		std_exons_log2_emb4 = np.std(exons_log2_emb4)

		gene_mean_std_exon_log2_emb4[gene] = [mean_exons_log2_emb4, std_exons_log2_emb4]


	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)
		
	U2_GTAG_5_win = range(13)
	U2_GTAG_5_max_score = 0
	index = 0
	for N in U2_GTAG_5_win:
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])
		index += 1
	
	U2_GTAG_3_win = range(17)
	U2_GTAG_3_max_score = 0
	index = 0
	for N in U2_GTAG_3_win:
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])
		index += 1



	for g in gene_introns.items():


		gene, introns = g

		#if gene in genes_intronic_EISA | genes_exonic_EISA:

		introns_log2_emb4 = set([])

		for intron in introns:

			intron_log2_emb4 = intronic_log2_emb4[intron]

			mean_exons_log2_emb4, std_exons_log2_emb4 = gene_mean_std_exon_log2_emb4[gene]

			exonic_z_score = z_score(intron_log2_emb4, mean_exons_log2_emb4, std_exons_log2_emb4)

			introns_log2_emb4.add((intron, intron_log2_emb4, exonic_z_score))


		for i in introns_log2_emb4:

				intron, log2_emb4, exonic_z_score = i

				try:

					ichrom, istart, iend, istrand = intron.split("_")

					istart = int(istart) - 10
					iend = int(iend) + 10 

					SJ_5 = Genome[ichrom][(istart-3):(istart+10)]
					SJ_3 = Genome[ichrom][(iend-14):(iend+3)]

					intron_seq = str(Genome[ichrom][istart:iend]).upper()

					SJ_5 = str(SJ_5).upper()
					SJ_3 = str(SJ_3).upper()
					
					if istrand == "-":
						SJ_5 = Genome[ichrom][(iend-10):(iend+3)]
						SJ_3 = Genome[ichrom][(istart-3):(istart+14)]

						SJ_5 = str(SJ_5.reverse_complement()).upper()
						SJ_3 = str(SJ_3.reverse_complement()).upper()

						intron_seq = str(Genome[ichrom][istart:iend].reverse_complement()).upper()


					GC_percent = percent( len(re.findall("[GC]", intron_seq)), len(intron_seq) )


					U2_GTAG_5_score = 0
					index = 0 
					for N in SJ_5:
						try:
							frec = U2_GTAG_5[N][index]
							U2_GTAG_5_score += frec
							index += 1
						except KeyError:
							pass

					U2_GTAG_3_score = 0
					index = 0 
					for N in SJ_3:
						try:
							frec = U2_GTAG_3[N][index]
							U2_GTAG_3_score += frec
							index += 1					 
						except KeyError:
							pass




					log2_emb4_other_introns = [x[1] for x in set(introns_log2_emb4) - set([(i)]) ]

					mean_log2_emb4_other_introns = np.mean(log2_emb4_other_introns)

					std_log2_emb4_other_introns = np.std(log2_emb4_other_introns)

					intronic_z_score = 	z_score(log2_emb4, mean_log2_emb4_other_introns, std_log2_emb4_other_introns)				

					U2_GTAG_total_score = percent( (U2_GTAG_5_score + U2_GTAG_3_score), (U2_GTAG_5_max_score + U2_GTAG_3_max_score) )

					U2_GTAG_5_score = percent( U2_GTAG_5_score, U2_GTAG_5_max_score )
					U2_GTAG_3_score = percent( U2_GTAG_3_score, U2_GTAG_3_max_score )

			
					# print intron, SJ_5, SJ_3, U2_GTAG_total_score


					if gene in genes_intronic_EISA: 

						print gene, len(introns_log2_emb4), intron, GC_percent, U2_GTAG_5_score, U2_GTAG_3_score,  U2_GTAG_total_score, intronic_z_score, exonic_z_score, "intronic_EISA", intron_CELE45.get(intron, "False"), (intron in DEXSeq_introns)

					elif gene in genes_exonic_EISA:

						print gene, len(introns_log2_emb4), intron, GC_percent, U2_GTAG_5_score, U2_GTAG_3_score, U2_GTAG_total_score,  intronic_z_score, exonic_z_score, "exonic_EISA", intron_CELE45.get(intron, "False"), (intron in DEXSeq_introns)

					else:

						print gene, len(introns_log2_emb4), intron, GC_percent, U2_GTAG_5_score, U2_GTAG_3_score, U2_GTAG_total_score, intronic_z_score, exonic_z_score, "null_EISA", intron_CELE45.get(intron, "False"), (intron in DEXSeq_introns)


				except ValueError:  #skip the split error of the piRNA sensor name 
					pass


if __name__ == '__main__':
	Genomictabulator(sys.argv[10])
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9] )

#python ~/my_src/PhD/Alper/intron_centric_analysis.py  PE100.introns.featureCounts.flat PE100.exons.featureCounts.flat SX2929_SX2930.EISA.txt ce10.RepeatMasker WS220_gurdon_curated_introns_genes_pEM975_filtered.gtf WS220_gurdon_curated_exons_genes_pEM975_filtered.gtf DEXSeq_old.out.old GT_AG_U2.donor GT_AG_U2.acceptor ../../../Genome/ce10/ce10.fa

#python ~/my_src/PhD/Alper/intron_centric_analysis.py  PE100.introns.featureCounts.flat PE100.exons.featureCounts.flat SX2929_SX2930.EISA.txt ce10.RepeatMasker WS220_gurdon_curated_introns_genes_pEM975_filtered.gtf WS220_gurdon_curated_exons_genes_pEM975_filtered.gtf	