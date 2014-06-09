import sys
from Bio import SeqIO

def Ndetector(fasta):
	
	f = open(fasta)
	for seq_record in SeqIO.parse(f, "fasta"):
		if "N" in str(seq_record.seq)[14:15]:
			print ">" + seq_record.id + '\n' + seq_record.seq 
	f.close()


if __name__ == '__main__':                                      
	Ndetector(sys.argv[1])

