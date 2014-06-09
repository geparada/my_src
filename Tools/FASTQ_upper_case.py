import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def main(fastq):
	for record in SeqIO.parse(fastq, "fastq"):
		Q = record.letter_annotations["phred_quality"]
		
		upperseq = SeqRecord( Seq(str(record.seq).upper()), id = record.id, description = "" )
		upperseq.letter_annotations["phred_quality"] = Q
		print upperseq.format("fastq")


if __name__ == '__main__':
	main(sys.argv[1])  
