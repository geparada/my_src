import sys
import csv
from Bio import SeqIO

def FF(fasta):
	
	f = open(fasta)
	for seq_record in SeqIO.parse(f, "fasta"):
		Ns = str(seq_record.seq).count('N')
		if 'N' in str(seq_record.seq)[14:16]:
#			pass           #ver los que se eliminan
#		elif Ns!=0:
			print seq_record.id, seq_record.seq, str(seq_record.seq)[:15].count('N'), str(seq_record.seq)[15:].count('N') ,Ns
	f.close()


if __name__ == '__main__':                                      
	FF(sys.argv[1])

