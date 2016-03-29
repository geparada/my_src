import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def main(EECC_fastq):

	f = open(EECC_fastq)

	for ERCC in SeqIO.parse(f, "fasta"):

		
		n=-1
		poliA_count=0

		while ERCC.seq[n]=="A":

			poliA_count+=1
			n-=1

		chrID = ERCC.id
		start = "1"
		end =  str(len(ERCC.seq)-poliA_count)
		strand = "+"
		gene_id = "gene_id " + '"' + ERCC.id + "_G" + '";'
		transcript_i = "gene_id " + '"' + ERCC.id + "_T" + '";'

		print "\t".join([chrID, "exon", start, end, ".", strand,  ".", gene_id, transcript_i])


if __name__ == '__main__':
	main(sys.argv[1])