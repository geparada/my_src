import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


def main (fastq):
	""" Toma el remanente del trimming para los paired end, y hace reverse complement a los reads que fueron secuenciados al reves """
	
	
	
	for record in SeqIO.parse(fastq, "fastq"):
		
		Q = record.letter_annotations["phred_quality"]

		if record.id[-2:]=="_1":
		
			upperseq = SeqRecord( record.seq.reverse_complement(), id = record.id, description = "" )
			upperseq.letter_annotations["phred_quality"] = Q[::-1]
			print upperseq.format("fastq"),
		
		else:
			upperseq = SeqRecord( record.seq, id = record.id, description = "" )
			upperseq.letter_annotations["phred_quality"] = Q			
			print upperseq.format("fastq"),	


if __name__ == '__main__':
	main(sys.argv[1])  	
