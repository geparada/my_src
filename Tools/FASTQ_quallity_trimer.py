import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def main(fastq,threshold):
	for record in SeqIO.parse(fastq, "fastq"):
		
		nt_error = 0
		L = len(record.seq)	
				
		while record.letter_annotations["phred_quality"][(-nt_error-1)] < threshold:
			nt_error += 1
			if nt_error==L:
				break
		if nt_error == 0:
			sequence = str(record.seq)
			Q = record.letter_annotations["phred_quality"]
		else:								
			sequence = str(record.seq)[:-nt_error]
			Q = record.letter_annotations["phred_quality"][:-nt_error]		
		
		upperseq = SeqRecord( Seq(sequence.upper()), id = record.id, description = "" )
		upperseq.letter_annotations["phred_quality"] = Q
		
		if len(sequence)>=50:
			print upperseq.format("fastq"),		


if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2]))  
