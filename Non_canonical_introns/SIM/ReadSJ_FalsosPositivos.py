import sys
import random 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def GenomicReadsGen(fasta, r1, r2, n):
	''' Simula reads de genomico
 	USO: python ReadSJ_FalsosPositivos.py ALLgenome.fa tamano_menor tamano_mayor ciclos '''


	SeqTable=[]
        f = open(fasta)

	print >> sys.stderr, "Loanding genome into RAM memory ..."

        for chrfa in SeqIO.parse(f, "fasta"):
                table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)	
			
	for x in range(n):

		for chr in SeqTable:
			startP = random.randrange(0, 250000000)                   
			longP = random.randrange(r1, r2+1)
			endP = startP + longP

			if len(chr[1])>endP:                                    #Para asegurarse de que el index no este fuera de rango
				seqP = chr[1][startP:endP].lower()
				nsP = seqP.count('n')
				if nsP==0:                                      #Para generar reads si NNNNNNNNNNNNNs
					print '>' + chr[0] + ':' + str(startP) + '+' + str(endP) + '\n' + seqP

	                startN = random.randrange(0, 250000000)
	                longN = random.randrange(r1, r2+1)
	                endN = startN + longN

	                if len(chr[1])>endN:
	                        seqN = chr[1][startN:endN][::-1].lower()
	                        nsN = seqN.count('n')
				if nsN==0:
					print '>' + chr[0] + ':' + str(startN) + '-' + str(endN) + '\n' + seqN.complement()

        f.close()


if __name__ == '__main__':
	GenomicReadsGen(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) )
