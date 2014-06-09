import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def main(fastq):
	
	for record in SeqIO.parse(fastq, "fastq"):
		seq = record.seq
		ID = record.id
		Q = record.letter_annotations["phred_quality"]		

		if str(record.seq[4:6]) =="CG":

			fastq_out = SeqRecord( seq, id = ID, description = "" )
			fastq_out.letter_annotations["phred_quality"] = Q
			
			print fastq_out.format("fastq"),			


if __name__ == '__main__':
	main(sys.argv[1])  
