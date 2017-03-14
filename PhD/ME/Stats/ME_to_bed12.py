import sys
import csv
import re




def main(ME_final):

	reader = csv.reader(open(ME_final), delimiter = '\t')

	header = reader.next()

	for row in reader:




		ME, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, min_P_MEs, total_ME, ME_P_value, Filtered = row
		

		ME_chrom, strand, ME_start, ME_end  = ME.split("_")

		ME_start= int(ME_start)
		ME_end = int(ME_end)

		

		for SJ in total_SJs.split(","):

			SJ_chr, SJ_start, SJ_end = re.findall(r"[\w']+", SJ)


			SJ_start= int(SJ_start)
			SJ_end = int(SJ_end)

			name = "_".join(map(str, ([len_micro_exon_seq_found,  sum_total_coverage])))


			print "\t".join(map(str, (ME_chrom, SJ_start-8, SJ_end+8, name, 0, strand, ME_start-8, ME_end+8, "0,0,0", 3, ",".join(["8", len_micro_exon_seq_found, "8" ]) , ",".join(map(str, [0, ME_start-SJ_start+8, SJ_end-SJ_start+8])) ) ) )




if __name__ == '__main__':
	main(sys.argv[1])