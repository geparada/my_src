import sys
import csv







def main(fasta_dusmasker, row_ME):

	read_info = {}



	for row in csv.reader(open(row_ME), delimiter = '\t'):

		read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found, micro_exons_coords = row


		read_info[read] = row





	ID = []

	for row in csv.reader(open(fasta_dusmasker), delimiter = ' '):

		if row[0][0] == ">":


			ID = row[0].strip(">")


		else:


			lc_start, minus, lc_end = row

			lc_start = int(lc_start)
			lc_end = int(lc_end)
			
			read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found, micro_exons_coords = read_info[ID]


			q_block_starts = map(int, q_block_starts.split(","))
			q_block_ends = map(int, q_block_ends.split(","))

			if len(q_block_starts)==2 and len(q_block_ends)==2:


				q_ME_start = q_block_ends[0]
				q_ME_end  = q_block_starts[-1]

				q_ME_seq = seq[q_ME_start:q_ME_end]


				relative_ME_lc_start = 0
				relative_ME_lc_end =  0


				if q_ME_start <= lc_start:

					relative_ME_lc_start =  lc_start - q_ME_start

				elif q_ME_start < lc_start < q_ME_end:

					relative_ME_lc_start = 0
				else:

					 relative_ME_lc_start =  lc_start - q_ME_end

 
				if q_ME_start <= lc_end:

					relative_ME_lc_end =  lc_end - q_ME_start

				elif q_ME_start < lc_end < q_ME_end:

					relative_ME_lc_end = 0
				else:

					 relative_ME_lc_end = lc_end - q_ME_end



				lc_seq = seq[int(lc_start) : int(lc_end) + 1]

				if DR_corrected_micro_exon_seq_found == "":     #There is an error at the DR_corrected_micro_exon_seq_found. Some of the values are "" (empty strings)
					DR_corrected_micro_exon_seq_found = micro_exon_seq_found


				DR_corrected_micro_exon_seq_found_comp = "".join(set(DR_corrected_micro_exon_seq_found)-set("N"))
				lc_seq_comp = "".join(set(lc_seq)-set("N"))

				if len(DR_corrected_micro_exon_seq_found_comp)>=2:
					DR_corrected_micro_exon_seq_found_comp="Other"

				if len(lc_seq_comp)>=2:
					lc_seq_comp="Other"

				if len(DR_corrected_micro_exon_seq_found_comp)>0 and len(lc_seq_comp)>0:


					print read, seq, DR_corrected_micro_exon_seq_found, DR_corrected_micro_exon_seq_found_comp, lc_start, lc_end, lc_seq, lc_seq_comp, relative_ME_lc_start, relative_ME_lc_end




			# print start, end, ID


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])