import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def main(fasta):
	for record in SeqIO.parse(fasta, "fasta"):
		seq = record.seq
		new_seq = str(seq[len(seq)/2:]) + str(seq[:len(seq)/2])
		
		print ">" + record.id
		print new_seq
		
if __name__ == '__main__':
	main(sys.argv[1]) 
