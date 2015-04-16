import sys
import csv
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SJ_Tags_seq = {}
SJ_Tags_info = {}

def Tagloader(fasta):
	
	print >> sys.stderr, "Loading SJ Tags in RAM memory ...",	

	f = open(fasta)

	for tag in SeqIO.parse(f, "fasta"):
		SJ_Tags_seq[tag.id.split("|")[0]] = tag.seq
		SJ_Tags_info[tag.id.split("|")[0]] = tag.id
		
	print >> sys.stderr, "OK"

	f.close()

def main(ME_centric_filter3):

	for row in csv.reader(open(ME_centric_filter3), delimiter = '\t'):

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_ME, total_ME, true_ME, score, is_annotated = row

		for SJ in total_SJs.split(","):

			SJ_Tag_seq = SJ_Tags_seq[SJ]
			up_block, down_block =  SJ_Tags_info[SJ].split("|")[-1].split("_")

			ME_Tag_ID = "|".join(SJ_Tags_info[SJ].split("|")[:-1] + ["_".join([up_block, len_micro_exon_seq_found, down_block])] )
			ME_Tag_seq = SJ_Tag_seq[:int(up_block)] + micro_exon_seq_found + SJ_Tag_seq[int(up_block):] 		

			print ">" + ME_Tag_ID
			print ME_Tag_seq

if __name__ == '__main__':
	Tagloader(sys.argv[1])
	main(sys.argv[2])	