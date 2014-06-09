import sys
from Bio import SeqIO

A=[]

def IDcompare(repeatmasker,dustmasker):
	
	f1 = open(repeatmasker)
	f2 = open(dustmasker)
	
	for RM, DM in zip(SeqIO.parse(f1, "fasta"), SeqIO.parse(f2, "fasta")):


		#if str(RM.seq.upper!=DM.seq.upper:
		#	print RM.id, RM.seq, DM.seq
		T = RM.id, RM.seq, DM.id, DM.seq
		print T			


#	fRM = SeqSeqIO.parse(f1, "fasta")
#	fDM = SeqIO.parse(f2, "fasta")
#	print fRM.id




	f1.close()
	f2.close()


if __name__ == '__main__':                                      
	IDcompare(sys.argv[1],sys.argv[2])

