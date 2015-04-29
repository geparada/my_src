import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
import re

#chr, start, end = re.findall(r"[\w']+", intron)

tags = {}


def Tags_indexer(tags_fasta):
	
	print >> sys.stderr, "Cargando a fasta en la ram ...",
	
	for record in SeqIO.parse(tags_fasta, "fasta"):
		if len(record.id.split("|")[-1].split("_")) == 3:
			id = record.id.split("|")[0]
			tags[id] = str(record.seq)
		
	print >> sys.stderr, "OK"


def main(ME_centric_filter3, gencode_bed12, ME_SJ_file):

	exon5_exon = defaultdict(set)
	exon3_exon = defaultdict(set)

	estart_introns = defaultdict(set)
	eend_introns = defaultdict(set)



	for row in csv.reader(open(gencode_bed12), delimiter = '\t'):
		
		csv.field_size_limit(1000000000)

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]


		#for q1, b in zip(qstarts, blocksizes):
		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):

			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = chr + ":" + str(istart) + strand + str(iend)


			estart = start + q1
			eend = start + q1 + b
			elenght = eend - estart
			exon = chr + ":" +  str(estart) + strand + str(eend)

			estart_introns[chr + "_" +  str(iend)].add(intron)
			eend_introns[chr + "_" +  str(istart)].add(intron)			

			if eend - estart > 25:



				if strand == "+":

					exon5_exon[chr + "_" +  str(estart)].add(exon)
					exon3_exon[chr + "_" +  str(eend)].add(exon)

				elif strand == "-":

					exon3_exon[chr + "_" +  str(estart)].add(exon)
					exon5_exon[chr + "_" +  str(eend)].add(exon)

	ME_SJ = defaultdict(int)
	
	for row in csv.reader(open(ME_SJ_file), delimiter = '\t'):

		read, flag, tag, start, cigar, seq, qual = row

		intron_tag, transcript_ID, anchors = tag.split("|")

		anchor_ME = 0
		ME_seq = ""

		if len(anchors.split("_"))==2:
			anchor_up, anchor_down = anchors.split("_")							

		if len(anchors.split("_"))==3:
			anchor_up, anchor_ME, anchor_down = anchors.split("_")


		anchor_up = int(anchor_up)
		anchor_ME = int(anchor_ME)
		anchor_down = int(anchor_down)

		if anchor_ME ==0:

			ME_SJ[intron_tag]+=1

		else:

			ME_seq = tags[intron_tag][anchor_up:anchor_up+anchor_ME]
			ME_SJ[intron_tag+"_"+str(anchor_ME)]+=1

		# if tag=="chr10:79797062+79800372|ENST00000372360.3|100_22_83":
		# 	print ME_seq


	for row in csv.reader(open(ME_centric_filter3), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_ME, total_ME, true_ME, score, is_annotated = row


		ME_SJ_coverages = []
		SJ_coverages = []
	

		ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = true_ME.split("|")
		ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")
		
		for SJ in total_SJs.split(","):
			SJ_coverage = ME_SJ[SJ]
			ME_SJ_coverage = ME_SJ[SJ+"_"+len_micro_exon_seq_found]

			#print SJ+"_"+micro_exon_seq_found

			ME_SJ_coverages.append(ME_SJ_coverage)
			SJ_coverages.append(SJ_coverage)

		sum_ME_coverage = sum(ME_SJ_coverages)
		sum_SJ_coverage = sum(SJ_coverages)

		ME_SJ_coverages = ",".join(map(str, ME_SJ_coverages))
		SJ_coverages = ",".join(map(str, SJ_coverages))


		is_alternative_5 = False
		is_alternative_3 = False

		alternatives_5 = []
		cov_alternatives_5 = []

		alternatives_3 = []
		cov_alternatives_3 = []

		total_cov_alternatives_5_pairs = set([])
		total_cov_alternatives_3_pairs = set([])


		if ME_strand == "+":

			if (ME_chrom + "_" + ME_start) in exon5_exon:

				is_alternative_3 = True

				for e in exon5_exon[ME_chrom + "_" + ME_start]:
					alternatives_3.append(e)

					chr, estart, eend = re.findall(r"[\w']+", e)
					cov_e = 0
					for i in estart_introns[chr + "_" +  str(estart)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_3_pairs.add((chr + "_" +  str(estart)+"|"+i, ME_SJ[i]))
					for i in eend_introns[chr + "_" +  str(eend)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_3_pairs.add((chr + "_" +  str(eend)+"|"+i, ME_SJ[i]))
					cov_alternatives_3.append(cov_e)

			if (ME_chrom + "_" + ME_end) in exon3_exon:

				is_alternative_5 = True

				for e in exon3_exon[ME_chrom + "_" + ME_end]:
					alternatives_5.append(e)

					chr, estart, eend = re.findall(r"[\w']+", e)
					cov_e = 0
					for i in estart_introns[chr + "_" +  str(estart)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_5_pairs.add((chr + "_" +  str(estart)+"|"+i, ME_SJ[i]))
					for i in eend_introns[chr + "_" +  str(eend)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_5_pairs.add((chr + "_" +  str(eend)+"|"+i, ME_SJ[i]))
					cov_alternatives_5.append(cov_e)


		elif ME_strand == "-":

			if (ME_chrom + "_" + ME_end) in exon5_exon:

				is_alternative_3 = True

				for e in exon5_exon[ME_chrom + "_" + ME_end]:
					alternatives_3.append(e)

					chr, estart, eend = re.findall(r"[\w']+", e)
					cov_e = 0
					for i in estart_introns[chr + "_" +  str(estart)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_3_pairs.add((chr + "_" +  str(estart)+"|"+i, ME_SJ[i]))
					for i in eend_introns[chr + "_" +  str(eend)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_3_pairs.add((chr + "_" +  str(eend)+"|"+i, ME_SJ[i]))
					cov_alternatives_3.append(cov_e)				

			if (ME_chrom + "_" + ME_start) in exon3_exon:

				is_alternative_5 = True

				for e in exon3_exon[ME_chrom + "_" + ME_start]:
					alternatives_5.append(e)

					chr, estart, eend = re.findall(r"[\w']+", e)
					cov_e = 0
					for i in estart_introns[chr + "_" +  str(estart)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_5_pairs.add((chr + "_" +  str(estart)+"|"+i, ME_SJ[i]))
					for i in eend_introns[chr + "_" +  str(eend)]:
						cov_e += ME_SJ[i]
						total_cov_alternatives_5_pairs.add((chr + "_" +  str(eend)+"|"+i, ME_SJ[i]))
					cov_alternatives_5.append(cov_e)

		total_cov_alternatives_5 = 0
		for p in total_cov_alternatives_5_pairs:
			total_cov_alternatives_5 += p[1]

		total_cov_alternatives_3 = 0
		for p in total_cov_alternatives_3_pairs:
			total_cov_alternatives_3 += p[1]

		if alternatives_5 == []:
			alternatives_5 = "None"
		else:
			alternatives_5 = ",".join(alternatives_5)	

		if cov_alternatives_5 == []:
			cov_alternatives_5 = "None"	
		else:
			cov_alternatives_5 = ",".join(map(str, cov_alternatives_5))	

		if alternatives_3 == []:
			alternatives_3 = "None"
		else:
			alternatives_3 = ",".join(alternatives_3)	

		if cov_alternatives_3 == []:
			cov_alternatives_3 = "None"	
		else:
			cov_alternatives_3 = ",".join(map(str, cov_alternatives_3))		

		#### Cambio de formato para paper ###


		# ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")

		# ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end
		# total_SJs = total_SJs.replace("+", "-")

		#if sum_ME_coverage>=3:

		print "\t".join(map(str, [ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, is_annotated, total_SJs, ME_SJ_coverages, sum_ME_coverage, SJ_coverages, sum_SJ_coverage, is_alternative_5, is_alternative_3, alternatives_5, cov_alternatives_5, total_cov_alternatives_5, alternatives_3, cov_alternatives_3,  total_cov_alternatives_3 ]))


#chr19:1105808+1106398_TGCCATCAAGTGGAACTTCACCAAG

#chr1_+_76256198_76256202



if __name__ == '__main__':
	Tags_indexer(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4] )

#python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round2/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3.ME_tags.fa /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 adipose1x75.sam.ME_SJ

#python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round2/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3.ME_tags.fa /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 adipose1x75.sam.pre_processed.filter1