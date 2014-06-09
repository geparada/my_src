import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def main(fasta):
	for record in SeqIO.parse(fasta, "fasta"):
		lft = [str(0), record.id, str(len(record.seq)), record.id, str(len(record.seq))]
		
		print "\t".join(lft)	
	
if __name__ == '__main__':
	main(sys.argv[1]) 
