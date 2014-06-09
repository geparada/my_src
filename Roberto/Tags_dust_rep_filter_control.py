import csv
import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord



def main(tags):
	
	for record in SeqIO.parse(tags, "fasta"):
		
		sequence = str(record.seq)
		
		sequence =  sequence[(len(sequence)/2):] + sequence[:(len(sequence)/2)]          #Para hacer el control

		
		upperseq = SeqRecord( Seq(sequence.upper()), id = record.id, description = "" )
		#upperseq.letter_annotations["phred_quality"] = Q

		print upperseq.format("fasta"), 		


if __name__ == '__main__':
	main(sys.argv[1])  
