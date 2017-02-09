import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	
     

def IntronExtractor(bed12):
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

	

	for row in csv.reader(open(bed12), delimiter = '\t'):
		
		csv.field_size_limit(1000000000)

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]


		transcript_length = sum(blocksizes)
		pos = 0 

		for q, b in zip(qstarts, blocksizes):
			estart = start + q
			eend = start + q + b
			elenght = eend - estart
			exon = chr + ":" +  str(estart) + strand + str(eend)


			for n in range(estart, eend):

				relative_pos = percent(pos,  transcript_length)

				if strand == "-":

					relative_pos = percent(transcript_length - pos,  transcript_length)

				pos += 1



				print "\t".join([exon, n, transcript_length])


				

			


if __name__ == '__main__':
	IntronExtractor(sys.argv[1])

