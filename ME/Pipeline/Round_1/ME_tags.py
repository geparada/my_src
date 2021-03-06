import sys
import csv
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from random import randint, sample
from operator import itemgetter
from collections import defaultdict
from operator import itemgetter


Transcriptome = {}

	
def Transcriptometabulator(genecode_fasta):
	
	print >> sys.stderr, "Cargando a fasta en la ram ...",
	
	for record in SeqIO.parse(genecode_fasta, "fasta"):
		id = str(record.id).split("|")[0]
		Transcriptome[id] = record.seq
		
	print >> sys.stderr, "OK"


def main(ME_centric_filter_3, bed12):

	n = 100

	transcript_intron_info = defaultdict(list)

	min_intron_lenght = 80

	SJ_ME = {}


	for row in csv.reader(open(ME_centric_filter_3), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, P_ME, total_ME, True_ME, True_ME_score, annotated  = row


		for i in total_SJs.split(","):  #Analizando las SJ host

			chrom, chromStart, chromEnd = re.findall(r"[\w']+", i)

			chromStart = int(chromStart)
			chromEnd = int(chromEnd)

			SJ_ME[i] = micro_exon_seq_found

			#print micro_exon_seq_found, chrom, chromStart, chromEnd



	for row in csv.reader(open(bed12), delimiter = '\t'):
		
		try:
		
			qName = row[3]
			seq = Transcriptome[qName]

			qstarts = map (int, row[11].strip(",").split(","))                      
			blocksizes = map(int, row[10].strip(",").split(","))

			start = int(row[1])
			strand = row[5]
			bn = int(row[9])
			chr = row[0]
			qstart = 0

			for q1, q2, b, b2 in zip(qstarts, qstarts[1:], blocksizes, blocksizes[1:]):
				
				qstart = qstart + b
				tag_start = qstart - n
				tag_end = qstart + n

				#if tag_start <= 0:
				#	print tag_start, qstart, tag_end, strand

				istart = start + q1 + b
				iend = start + q2
				ilen = iend - istart
				intron = row[0] + ":" +  str(istart) + row[5] + str(iend)	
				intron = chr + ":" + str(istart) + strand + str(iend)
				ilength = iend - istart

				block_up = n
				block_down = n
				
				if strand == '+' :                          #Para los que aliniean en la hebra +
								   
					if tag_start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
						tag_start = 0
						block_up = qstart

					if tag_end>len(seq):
						tag_end=len(seq)
						block_down = tag_end - qstart


					tag = seq[tag_start:tag_end]
					
								  
				if strand == '-' :
				
					if tag_end>len(seq):                 #Para los que alinian en la hebra - es todo al inverso
						tag_end=len(seq)
						block_up = tag_end - qstart

					tag = seq[-tag_end:-tag_start]

					if tag_start<=0:

						tag = seq[-tag_end:]
						block_down = qstart

										 
				if b > 25 and b2 > 25 and ilength >= min_intron_lenght:  # hay que agregarle el filtro de los micro exones!!

					info = qName, tag, chr, istart, iend, strand, block_up, block_down, block_up + block_down
					transcript_intron_info[intron].append(info)	


		except KeyError:
			pass


	for i in transcript_intron_info.items():

		infos = i[1]
		intron = i[0]

		qName, tag, chr, istart, iend, strand, block_up, block_down, sum_blocks = max(infos, key=itemgetter(8))


		ID = ">" + intron + "|" + qName + "|" + str(block_up) + "_" + str(block_down)

		#print ID
		#print tag

		if intron in SJ_ME:


			ME = SJ_ME[intron]
			ME_len = len(ME)

			new_ID = ">" + intron + "|" + qName + "|" + str(block_up) + "_" + str(block_down) + "|" + ME
			new_tag = tag[:block_up] + ME +  tag[block_up:]

			print new_ID
			print new_tag


if __name__ == '__main__':
	Transcriptometabulator(sys.argv[1])
	main (sys.argv[2],sys.argv[3]) 		


#python ~/my_src/ME/Pipeline/ME_tags.py ~/db/transcriptome/hg19/Gene_models/gencode/v11/gencode.v11.pc_transcripts.fa  TOTAL.filter1.ME_centric.filter2.filter3 ~/db/transcriptome/hg19/Gene_models/gencode/v11/gencode.v11.annotation.bed12