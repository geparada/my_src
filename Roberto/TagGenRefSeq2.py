import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



RNA = []
Genoma = []
SeqName = []
Resultados=[]


Genome = {}
Trascriptome = {}


def Genomictabulator(fasta):

	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	print >> sys.stderr, "OK"

	f.close()
	
def Transcriptometabulator(refseq_fasta):
	
	print >> sys.stderr, "Cargando transcriptoma en la memoria RAM ...",	
	
	for record in SeqIO.parse(refseq_fasta, "fasta"):
		Trascriptome[record.id] = record.seq
		
	print >> sys.stderr, "OK"	

def main(refseq_psl,len_tag):
	'''Generates tags of the splice-juntions
	*Input file = USCS BLAT alignments table (22 columns ) joined with fasta converted to tab (3 columns). This file has 25 columns.
	*Usage = python Tag-Fulinfo inputfile tag length
	*Output = seq ID , intron name, strand, tag sequence, tag length '''
       
   
	for row in csv.reader(open(refseq_psl), dialect='excel-tab' ):
		
		qstarts = map (int, row[20].split(",")[1:-1])                      # obteniendo lista int de las distintas variables
		tstarts = map(int, row[21].split(",")[0:-1])
		blocksizes = map(int, row[19].split(",")[0:-1])

#		print row, blocksizes, qstarts, tstarts 

		n = len_tag/2                                                      #preparando argumento para su uso

        
		for t, b, t2, b2, q in zip(tstarts, blocksizes, tstarts[1:], blocksizes[1:]  , qstarts):
			tbsum = t+b                                      # generando primera coordenada del intron al cual se le va a hacer el tag              
			start = q-n                                     #generando intervalos para hacer el slicing del la secuencia                      
			end = q+n
			chr = row[14]
			strand = row[9]
			qName = row[10]
			
#			print row
			
                    
			if strand == '+' :                          #Para los que aliniean en la hebra +
    
                    
				if start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
					start = 0  
				if end>len(Trascriptome[qName]):
					end=len(Trascriptome[qName])
				tag = Trascriptome[qName][start:end]
    
                    
			if strand == '-' :                                #Para los que alinian en la hebra - es todo al inverso
    
				if end>len(Trascriptome[qName]):
					end=len(Trascriptome[qName])                                       
				tag = Trascriptome[qName][-end:-start]
				if start<=0: 
					tag = Trascriptome[qName][-end:]
                         
			istart = tbsum
			iend = t2
			ilength = iend - istart
				
			intron = chr + ":" + str(istart) + strand + str(iend)
			refseq = row[10]
			
			dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]


			if strand == '-':
					dn = dn.reverse_complement()
					
			dn = str(dn)

					
			if (dn=="GTAG" or dn=="GCAG" or dn=="ATAC") and ilength >= 40 and len(tag)==len_tag:
					print ">" + refseq + "|" + intron + "\n" + tag
				
			
				


    
    
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	Transcriptometabulator(sys.argv[2])                
	main(sys.argv[3], int(sys.argv[4]))


#SJtag("BLAToincDNAmm9",  30, "BLATIntronescDNAmm9")
#switch("BLATIntronescDNAmm9")
#CompareTags("BLATIntronescDNAmm9")
