
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna




Genoma = []
SeqTable = []


	
def Transcriptometabulator(mrna_fasta):

	print >> sys.stderr, "Cargando fasta en la ram ...",
		
	for record in SeqIO.parse(mrna_fasta, "fasta"):
		table = (record.id, record.seq)
		SeqTable.append(table)
		
	print >> sys.stderr, "OK"

		
	
def main(psl, len_tag):
	'''Generates tags of the splice-juntions
	*Input file = USCS BLAT alignments table (22 columns ) joined with fasta converted to tab (3 columns). This file has 25 columns.
	*Usage = python Tag-Fulinfo inputfile tag length
	*Output = seq ID , intron name, strand, tag sequence, tag length '''
       
	print >> sys.stderr, "Extrayendo intrones del psl ...",
    
	reader1 = csv.reader(open(psl), dialect='excel-tab' )

         
         # row[0]  matches
         # row[1]  misMatches
         # row[2]  repMaches
         # row[3]  nCount
         # row[4]  qNumInsert    
         # row[5]  qBaseInsert
         # row[6]  tNumInsert
         # row[7]  tBaseInsert
         # row[8]  strand
         # row[9] qName
         # row[10] qSize
         # row[11] qStart
         # row[12] qEnd         
         # row[13] tName
         # row[14] tSize
         # row[15] tStart
         # row[16] tEnd
         # row[17] blockCount
         # row[18] blockSizes
         # row[19] qStarts
         # row[20] tStarts


   
	for row in reader1:

		try:
		
			Transcriptome = dict(SeqTable)
			qName = row[9]
			seq = Transcriptome[qName]
					
			qstarts = map (int, row[19].split(",")[1:-1])                      # obteniendo lista int de las distintas variables
			tstarts = map(int, row[20].split(",")[0:-1])
			blocksizes = map(int, row[18].split(",")[0:-1])
		

			n = len_tag/2                                                      #preparando argumento para su uso

			
			for t, b, t2, b2, q in zip(tstarts, blocksizes, tstarts[1:], blocksizes[1:]  , qstarts):
				tbsum = t+b                                      # generando primera coordenada del intron al cual se le va a hacer el tag              
				start = q-n                                     #generando intervalos para hacer el slicing del la secuencia                      
				end = q+n
				chr = row[13]
				strand = row[8]
				

				istart = tbsum
				iend = t2
				ilength = iend - istart				
				intron = chr + ":" + str(istart) + strand + str(iend)
								
				if strand == '+' :                          #Para los que aliniean en la hebra +
								   
					if start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
						start = 0  
					if end>len(seq):
						end=len(seq)
					tag = seq[start:end]
								  
				if strand == '-' :                                #Para los que alinian en la hebra - es todo al inverso
				
					if end>len(seq):
						end=len(seq)                                       
					tag = seq[-end:-start]
					if start<=0: 
						tag = seq[-end:]
										 
				if ilength >= 40:
						print intron, tag
						
		except KeyError:
			pass

	print >> sys.stderr, "OK"
			
					
			  
if __name__ == '__main__':
	Transcriptometabulator(sys.argv[1])
	main(sys.argv[2], int(sys.argv[3]))
	          



