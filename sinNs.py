import sys
import csv
from Bio import SeqIO

def FF(fasta):
	
	f = open(fasta)
	for seq_record in SeqIO.parse(f, "fasta"):
		Ns = str(seq_record.seq).count('N')
		if Ns==0:
			print ">" + seq_record.id + '\n' + seq_record.seq 
	f.close()


if __name__ == '__main__':                                      
	FF(sys.argv[1])

