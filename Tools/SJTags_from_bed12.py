import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna




Genoma = []
SeqTable = []


	
def Transcriptometabulator(genecode_fasta):
	
	print >> sys.stderr, "Cargando a fasta en la ram ...",
	
	for record in SeqIO.parse(genecode_fasta, "fasta"):
		id = str(record.id).split("|")[0]
		table = (str(id), record.seq)
		SeqTable.append(table)
		
	print >> sys.stderr, "OK"
	
		
	
def main(bed12, len_tag):

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

	Transcriptome = dict(SeqTable)


	for row in reader1:
		
		try:
		
			qName = row[3]
			seq = Transcriptome[qName]

			qstarts = map (int, row[11].strip(",").split(","))                      
			blocksizes = map(int, row[10].strip(",").split(","))

			start = int(row[1])
			strand = row[5]
			bn = int(row[9])
			chr = row[0]
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
				
				if strand == '+' :                          #Para los que aliniean en la hebra +
								   
					if tag_start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
						tag_start = 0  
					if tag_end>len(seq):
						tag_end=len(seq)
					tag = seq[tag_start:tag_end]
								  
				if strand == '-' :                                #Para los que alinian en la hebra - es todo al inverso
				
					if tag_end>len(seq):
						tag_end=len(seq)                                       
					tag = seq[-tag_end:-tag_start]
					if tag_start<=0: 
						tag = seq[-tag_end:]
										 
				if ilength >= 40:
						print intron, tag
						
		except KeyError:
			pass

	print >> sys.stderr, "OK"

			
			  
if __name__ == '__main__':
	Transcriptometabulator(sys.argv[1])
	main(sys.argv[2], int(sys.argv[3]))
             

