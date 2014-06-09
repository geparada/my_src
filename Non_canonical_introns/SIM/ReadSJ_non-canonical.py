import sys
import random 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def GenomicReadsGen(fasta, r1, r2, n):
	SeqTable=[]
        f = open(fasta)

        for chrfa in SeqIO.parse(f, "fasta"):
                table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)	
			
	for x in range(n):

		for chr in SeqTable:
			startP1 = random.randrange(0, 250000000)
			longP1 = random.randrange(1, 100)
			endP1 = startP1 + longP1
			
			intron_lengthP = random.randrange(40, 10000)

			startP2 = endP1 + intron_lengthP			
			longP2 = random.randrange(1, 100)
			endP2 = startP2 + longP2

			if len(chr[1])>endP2 and (longP1+longP2)>=r1 and (longP1+longP2)<=r2:
				seqP = chr[1][startP1:endP1].lower() + chr[1][startP2:endP2].lower()
				nsP = seqP.count('n')
				if nsP==0: 
					print '>' + chr[0] + ':' + str(endP1) + '+' + str(startP2) + '=' + str(longP1) + '-' + str(longP2) + '\n' + seqP         
 
			startN1 = random.randrange(0, 250000000)
			longN1 = random.randrange(1, 100)
			endN1 = startN1 + longN1
			
			intron_lengthN = random.randrange(40, 10000)

			startN2 = endN1 + intron_lengthN			
			longN2 = random.randrange(1, 100)
			endN2 = startN2 + longN2

			if len(chr[1])>endN2 and (longN1+longN2)>=r1 and (longN1+longN2)<=r2:
				seqN = chr[1][startN1:endN1].lower() + chr[1][startN2:endN2].lower()
				nsN = seqN.count('n')
				if nsN==0:
					print '>' + chr[0] + ':' + str(endN1) + '-' + str(startN2) + '=' + str(longN1) + '-' + str(longN2) + '\n' + seqN.reverse_complement()  




#	                if len(chr[1])>endN:
#	                        seqN = chr[1][startN:endN][::-1].lower()
#	                        nsN = seqN.count('n')
#				if nsN==0:
#					print '>' + chr[0] + ':' + str(startN) + '-' + str(endN) + '\n' + seqN.complement()

        f.close()


if __name__ == '__main__':
	GenomicReadsGen(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) )
