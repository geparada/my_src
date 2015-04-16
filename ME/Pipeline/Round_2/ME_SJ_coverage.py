import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict

tags = {}


def Tags_indexer(tags_fasta):
	
	print >> sys.stderr, "Cargando a fasta en la ram ...",
	
	for record in SeqIO.parse(tags_fasta, "fasta"):
		if len(record.id.split("|")[-1].split("_")) == 3:
			id = record.id.split("|")[0]
			tags[id] = str(record.seq)
		
	print >> sys.stderr, "OK"


def main(ME_centric_filter3, ME_SJ_file):

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
		total_number_of_micro_exons_matches

		ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates = true_ME.split("|")
		
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

		#### Cambio de formato para paper ###


		# ME_chrom, ME_strand, ME_start, ME_end = ME.split("_")

		# ME_ID = ME_chrom + ":" + ME_start + "-" + ME_end
		# total_SJs = total_SJs.replace("+", "-")

		#if sum_ME_coverage>=3:

		print "\t".join(map(str, [ME, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, is_annotated, total_SJs, ME_SJ_coverages, sum_ME_coverage, SJ_coverages, sum_SJ_coverage ]))


#chr19:1105808+1106398_TGCCATCAAGTGGAACTTCACCAAG





if __name__ == '__main__':
	Tags_indexer(sys.argv[1])
	main(sys.argv[2], sys.argv[3] )

#python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round2/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3.ME_tags.fa /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 adipose1x75.sam.ME_SJ