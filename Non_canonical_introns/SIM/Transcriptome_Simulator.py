import sys
import random 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import csv

         # row[0]  matches
         # row[1]  misMatches
         # row[2]  repMaches
         # row[3]  nCount
         # row[4]  qNumInsert    
         # row[5]  qBaseInsert
         # row[6]  tNumInsert
         # row[7]  tBaseInsert
         # row[8]  strand   O <------ 
         # row[9] qName     O <------
         # row[10] qSize
         # row[11] qStart
         # row[12] qEnd         
         # row[13] tName   O <------
         # row[14] tSize   
         # row[15] tStart  I <-----
         # row[16] tEnd    I <-----
         # row[17] blockCount  I,O <-----
         # row[18] blockSizes  O <-----
         # row[19] qStarts  O <----
         # row[20] tStarts  O <----

         # row[21] qName O <-----
         # row[22] 1
         # row[23] sequence  O <-----



def TranscriptomeGen(fasta, psl, n):
"""Randomiza un transcriptoma input """

	SeqTable=[]
        f = open(fasta)
	
	reader = csv.reader(open(psl), dialect='excel-tab' )


        for chrfa in SeqIO.parse(f, "fasta"):
                table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)


        for cycle in range(n):

		for row in reader:	
			

			t1 = [int(row[15])]
			tz = int(row[16])
			bc = int(row[17])
			
			b = [random.randrange(10,300)]
			q = [0]
						

			for i in range(320):                      #numero maximo de bloques (en el trascriptoma real de hg19 es 312)
				t2min = t1[-1] + 390              #El exon mas grande prodrá ser de 350 pb
				t2max = t1[-1] + 10000
				t2 = random.randrange(t2min,t2max)

				if t2<=tz :
					t1 = t1 + [t2]
				else:
					break				
				

			for i in range(len(t1)-1):
				b = b + [random.randrange(10,351)]
				q = q + [b[-2]+q[-1]]


                        tstarts = ''
			blocksizes = ''
                        qstarts = ''

			for i in range(len(t1)):
		
				tstarts = tstarts + str(t1[i]) + ','
				blocksizes = blocksizes + str(b[i]) + ','
				qstarts = qstarts + str(q[i]) + ','
				blocknum = str(len(t1))

			for chr in SeqTable:

				if chr[0] == row[13]:
					refseq = ''					
					for i in range(len(t1)):
						start = t1[i]
						end = t1[i] + b[i]	
						refseq = refseq + chr[1][start:end]
			#		print ">" + row[9] + "\n" + refseq 
										
 
			if row[8] == '+' :
											
				restult = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7],  row[8], row[9], row[10], row[11], str(q[-1]+b[-1]), row[13], row[14], row[15], str(int(t1[-1])+int(b[-1])), blocknum,  blocksizes,  qstarts, tstarts, row[9], 1, refseq.lower() ]

			if row[8] == '-' :

				restult = [ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7],  row[8], row[9], row[10], row[11], str(q[-1]+b[-1]), row[13], row[14], row[15], str(int(t1[-1])+int(b[-1])), blocknum,  blocksizes,  qstarts, tstarts, row[9], 1, refseq.lower().reverse_complement() ]
							 

			c = map(str, restult)
			
			print c[0] + '\t' + c[1] + '\t' + c[2] + '\t' + c[3] + '\t' + c[4] + '\t' + c[5] + '\t' + c[6] + '\t' + c[7] + '\t' + c[8] + '\t' + c[9] + '\t' + c[10] + '\t' + c[11] + '\t' + c[12] + '\t' + c[13] + '\t' + c[14] + '\t' + c[15] + '\t' + c[16] + '\t' + c[17] + '\t' + c[18] + '\t' + c[19] + '\t' + c[20] + '\t' + c[21] + '\t' + c[22] + '\t' + c[23]
			








		f.close()




if __name__ == '__main__':
	TranscriptomeGen(sys.argv[1], sys.argv[2], int(sys.argv[3]) )
