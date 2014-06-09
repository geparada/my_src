import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna




Genoma = []
SeqTable = []


	
def Transcriptometabulator(genecode_fasta, mrna_fasta, EST_fasta):
	
	print >> sys.stderr, "Cargando a GENCODE en la ram ...",
	
	for record in SeqIO.parse(genecode_fasta, "fasta"):
		id = str(record.id).split("|")[0]
		table = (str(id), record.seq)
		SeqTable.append(table)
		
	print >> sys.stderr, "OK"
	
	print >> sys.stderr, "Cargando mRNAs en la ram ...",
		
	for record in SeqIO.parse(mrna_fasta, "fasta"):
		table = (record.id, record.seq)
		SeqTable.append(table)
		
	print >> sys.stderr, "OK"
	print >> sys.stderr, "Cargando EST en la ram ...",	
		
	for record in SeqIO.parse(EST_fasta, "fasta"):
		table = (record.id, record.seq)
		SeqTable.append(table)
	
	print >> sys.stderr, "OK"
		
	

def IntronExtractor_bed12(bed12, len_tag):

	print >> sys.stderr, "Extrayendo intrones del bed12 ...",

         # row[0]  chr
         # row[1]  alignment start
         # row[2]  alignment end
         # row[3]  name 
         # row[4]    
         # row[5]  strand
         # row[6]  aligment start
         # row[7]  aligment end
         # row[8]  
         # row[9] blocknum
         # row[10] blocksizes
         # row[11] qstarts  


	
	reader1 = csv.reader(open(bed12), delimiter = '\t')

	csv.field_size_limit(1000000000)	
	
	n = len_tag/2 

	for row in reader1:

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]
		qName = row[3]
		qstart = 0

		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
			
			qstart = qstart + b
			tag_start = qstart - n
			tag_end = qstart + n

			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = row[0] + ":" +  str(istart) + row[5] + str(iend)	
			intron = chr + ":" + str(istart) + strand + str(iend)
			ilength = iend - istart

			Trascriptome = dict(SeqTable)
					
			try:
				
				if strand == '+' :                          #Para los que aliniean en la hebra +
								   
					if tag_start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
						tag_start = 0  
					if tag_end>len(Trascriptome[qName]):
						tag_end=len(Trascriptome[qName])
					tag = Trascriptome[qName][tag_start:tag_end]
								  
				if strand == '-' :                                #Para los que alinian en la hebra - es todo al inverso
				
					if tag_end>len(Trascriptome[qName]):
						tag_end=len(Trascriptome[qName])                                       
					tag = Trascriptome[qName][-tag_end:-tag_start]
					if tag_start<=0: 
						tag = Trascriptome[qName][-tag_end:]
										 
				if ilength >= 40:
						print intron, tag
						
			except KeyError:
				pass

			
	print >> sys.stderr, "OK"
	

def IntronExtractor_psl(psl, len_tag):
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
			qName = row[9]

			istart = tbsum
			iend = t2
			ilength = iend - istart				
			intron = chr + ":" + str(istart) + strand + str(iend)
			
			Trascriptome = dict(SeqTable)
			
			try:
				
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
										 
				if ilength >= 40:
						print intron, tag
						
			except KeyError:
				pass

	print >> sys.stderr, "OK"
			
					
			  
if __name__ == '__main__':
	Transcriptometabulator(sys.argv[1],sys.argv[2],sys.argv[3])
	IntronExtractor_bed12(sys.argv[4], int(sys.argv[7]))
	IntronExtractor_psl(sys.argv[5], int(sys.argv[7]))
	IntronExtractor_psl(sys.argv[6], int(sys.argv[7]))                



