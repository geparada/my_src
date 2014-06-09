import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord



 

def main(tags_sam):
	"""Filtra aquellos reads que alinearon a mas de un SJ """
	
	reader1 = csv.reader(open(tags_sam), delimiter = '\t')
	reader2 = csv.reader(open(tags_sam), delimiter = '\t')	
	
	intron_reads = defaultdict(set)
	output_reads = set([])
	
	for row in reader1:
		
		if row[1]=='0' or row[1]=='16':
	
			read = row[0]		
			flag = int(row[1])
			tag_name = row[2]
			start = int(row[3])           #Sam es 1 referenciado 
			cigar = row[5]
			seq = row[9]
			qual = row[10]
			intron_tag = tag_name.split("|")[0]				
		
			intron_reads[read].add(intron_tag)

	for row in reader2:
		
		if row[1]=='0' or row[1]=='16':
		
			read = row[0]		
			flag = int(row[1])
			tag_name = row[2]
			start = int(row[3])           #Sam es 1 referenciado 
			cigar = row[5]
			seq = row[9]
			qual = row[10]
			intron_tag = tag_name.split("|")[0]
			
			if len(intron_reads[read]) == 1:
				fastq = ("@" + read + "|" + intron_tag, seq, "+", qual)
				output_reads.add(fastq)

	for r in output_reads:
		print r[0]
		print r[1]
		print r[2]
		print r[3]

#			fastq = SeqRecord( Seq(seq), id = read + "|" + intron_tag, description = "" )
#			fastq.letter_annotations["phred_quality"] = qual
#			print fastq.format("fastq")
			

	

if __name__ == '__main__':
	main(sys.argv[1])
