import sys
from Bio import SeqIO

def IDfinder(fasta,ID):
	
	f = open(fasta)
	for seq_record in SeqIO.parse(f, "fasta"):
		if seq_record.id==ID:
			print ">" + seq_record.id + '\n' + seq_record.seq 
	f.close()


if __name__ == '__main__':                                      
	IDfinder(sys.argv[1],sys.argv[2])

