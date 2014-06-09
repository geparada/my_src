import sys
import csv
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#

RNA = []
Genoma = []
SeqName = []
Resultados=[]
Transcriptome = []


SeqTable = []


def Genomictabulator(fasta):
	
	f = open(fasta)


	for chrfa in SeqIO.parse(f, "fasta"):
		table = (str(chrfa.id), chrfa.seq)
		SeqTable.append(table)

	f.close()


def expression (gene_model):

	Genome = dict(SeqTable)

	reader1 = csv.reader(open(gene_model), delimiter='\t' )
		 
		# row[0]  chr
		# row[1]  tstart
		# row[2]  tend
		# row[3]  name
		# row[4]    
		# row[5]  strand
		# row[6]  tstart
		# row[7]  tend
		# row[8]  
		# row[9] blocknum
		# row[10] blocksizes
		# row[11] qstarts         
		
		# row[12] Seq
 

	for row in reader1:
		
		blocknum = int(row[9])
		gene = [row] + ['seq']
		#print gene

		if blocknum>=4:
			Transcriptome.append(gene)
	


def SJtag(phred_quality):
	Genome = dict(SeqTable)

	'''Genera reads con splice juntios desde un transcriptoma de referencia

	*Input file = USCS BLAT alignments table (21 columns ) joined with fasta converted to tab (3 columns). This file has 24 columns.
	*Usage = python file longitud_minima longitud_maxima numero_de_ciclos
	*Output = seq ID , intron name, strand, tag sequence, tag length '''

	reader2 = csv.reader(open(phred_quality), delimiter=' ' )
	
	for quality in reader2:
	

		gene = random.choice(Transcriptome)
		row, seq = gene

		start = int(row[1])

		qstarts = map (int, row[11].strip(",").split(","))
		tstarts = [x + start for x in qstarts]
		blocksizes = map(int, row[10].strip(",").split(","))


		strand = row[5]
		chromosome = row[0]
		gene_name = row[3]	
	       

		x = random.randint(0,len(qstarts)-4)
	
		t0 = tstarts[x]
		b0 = blocksizes[x]

		t1 = tstarts[x+1]
		b1 = blocksizes[x+1]

		t2 = tstarts[x+2]
		b2 = blocksizes[x+2]

		t3 = tstarts[x+3]
		b3 = blocksizes[x+3]

		q = qstarts[x+2]	


		l = len(quality[0])

			
		n1 = random.randint(1,l-1)
		n2 = l - n1		
			
	    
		tbsum0 = t0+b0
		tbsum1 = t1+b1  
		tbsum2 = t2+b2
			#tbsum3 no es necesario 
					                                           
		istart = q-n1                             
		iend = q+n2
		if istart>0 or iend>len(seq) or (b0+b1)<n1 or (b2+b3)<n2:     
			pass

	# Y que el tag no se pase de la zona del transcrito asegura que solo tengan 3 SJ como maximo


		intron0 = ''                                                  
		intron1 = chromosome + ":" + str(tbsum1) + strand + str(t2)                  
		intron2 = ''    

		anchor0 = ''
		anchor1 = str(n1) + "-" + str(n2)
		anchor2 = ''
		
		len0 = t1-tbsum0 
		len1 = t2-tbsum1
		len2 = t3-tbsum2 

		tag = Genome[chromosome][tbsum1-n1:tbsum1] + Genome[chromosome][t2:t2+n2]

		

	# los intrones 0 y 2 estaran presente solo si el start o el end se escapan del rango en el cual se encuentran los exones mas proximos al SJ

		if n1>b1 and n2<=b2:
			intron0 = chromosome + ":" + str(tbsum0) + strand + str(t1) + "<>"


			anchor0 = str(n1-b1) + "-" + str(b1) + "<>"
			anchor1 = str(b1) + "-" + str(n2)

			tag = Genome[chromosome][tbsum0-(n1-b1):tbsum0] + Genome[chromosome][t1:tbsum1] + Genome[chromosome][t2:t2+n2]

		if n2>b2 and n1<=b1:
			intron2 = "<>" + chromosome + ":" + str(tbsum2) + strand + str(t3)



			anchor1 = str(n1) + "-" + str(b2)
			anchor2 = "<>" + str(b2) + "-" + str(n2-b2)
			tag = Genome[chromosome][tbsum1-n1:tbsum1] + Genome[chromosome][t2:tbsum2] + Genome[chromosome][t3:t3+(n2-b2)]


		if n1>b1 and n2>b2:
			intron0 = chromosome + ":" + str(tbsum0) + strand + str(t1) + "<>"
			intron2 = "<>" + chromosome + ":" + str(tbsum2) + strand + str(t3)

			anchor0 = str(n1-b1) + "-" + str(b1) + "<>"
			anchor1 = str(b1) + "-" + str(b2)
			anchor2 = "<>" + str(b2) + "-" + str(n2-b2)
			tag = Genome[chromosome][tbsum0-(n1-b1):tbsum0] + Genome[chromosome][t1:tbsum1] + Genome[chromosome][t2:tbsum2] + Genome[chromosome][t3:t3+(n2-b2)]


		
		read = tag.lower()
            	            
		if strand == '-' :                            
			read = str(read.reverse_complement())


		print "@" + gene_name + "=" + intron0 + intron1 + intron2 + "=" + anchor0 + anchor1 + anchor2 + "\n" + read + "\n" + "+" + "\n" + quality[0]




    
if __name__ == '__main__':
	#print 'Cargando genoma en la ram...'
	Genomictabulator(sys.argv[1])
	#print 'Cargando modelo de transcriptoma en la ram ...'
	expression(sys.argv[2])                       
	SJtag(sys.argv[3])
