import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable = []

def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

	f.close()

     

def compare(clip_trim, tablasfinales, columna_intrones):
         # row[0]  chr
         # row[1]  alignment start
         # row[2]  alignment end
         # row[3]  name (Refseq=intron0<>intron1=anchor0<>anchor1)
         # row[4]    
         # row[5]  strand
         # row[6]  aligment start
         # row[7]  aligment end
         # row[8]  
         # row[9] blocknum
         # row[10] blocksizes
         # row[11] qstarts  


	reader1 = csv.reader(open(clip_trim), delimiter = '\t')
	reader2 = csv.reader(open(tablasfinales), delimiter = ' ')  #Depende si se usa genecode o tablas finales
	csv.field_size_limit(1000000000)	

	known = float(0)
	unknown = float(0)
	

	introns_known = []

	for row in reader2:
		introns_known.append(row[columna_intrones-1])	

	for row in reader1:
        	Genome = dict(SeqTable)

		#introns_set = set(introns_known)
		qstarts = map (int, row[11].strip(",").split(","))                      # obteniendo lista int de las distintas variables de la linea psl
		blocksizes = map(int, row[10]..strip(",").split(","))
		start = int(row[1])
		strand = row[5]
		bn = int(row[9])

		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = row[0] +  str(istart) + row[5] + str(iend)
			chr = row[0]


			result = ''

			if intron in introns_known:
				result = 'Known'
			else:
				result = 'Unknown'			
			

			dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
			if row[5] == '-':
					dn = dn.reverse_complement()
			if ilen >= 40 and result == 'Known' and str(dn).upper() != 'GTAG' and str(dn).upper() != 'GCAG' and str(dn).upper() != 'ATAC' :	
				print row[3], intron, result, str(dn).upper()			
			

			
		
					
				
			#print row[3], intron, result, dn

	#total = known + unknown	
	#print 'total =', total, 'Known = ', (100*known)/total, 'Unknown =', (100*unknown)/total
	  		

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	compare(sys.argv[2], sys.argv[3], int(sys.argv[4]))

