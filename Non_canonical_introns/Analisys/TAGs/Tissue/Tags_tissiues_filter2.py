import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main (fastq, repbase, dust_IDs, bowtie_genome):
	""" Filtra aquellos reads que pueden alinear al genoma, zonas repetitivas o son de baja complejidad """
	
	repeat_low_complex_genome = set([])
	
	intron_coverage = defaultdict(set)
		

	reader2 = csv.reader(open(repbase), delimiter =  '\t')
	reader3 = csv.reader(open(dust_IDs), delimiter =  ' ')
	reader4 = csv.reader(open(bowtie_genome), delimiter =  '\t')	
				
	
	for row in reader2:
				
		matches = int(row[0])
		qName = row[9]
		qSize = int(row[10])
		rep_percent = percent(matches, qSize)
		
		if rep_percent > 50:
			repeat_low_complex_genome.add(qName)
			
	for row in reader3:
		ID = row[0][1:]
		repeat_low_complex_genome.add(ID)
	
	for row in reader4:
		if row[1]=="0" or row[1]=="16":
			repeat_low_complex_genome.add(row[0])				
			
				
	for record in SeqIO.parse(fastq, "fastq"):
		read = record.id.split("|")[0]
		intron = record.id.split("|")[1]
		seq = str(record.seq) 
		
		if not record.id in repeat_low_complex_genome:
			intron_coverage[intron].add(read + "|" + seq )
	
	for i in intron_coverage.items():
		introns = i[0]
		#reads = list(i[1])
		coverage = len(i[1])
		reads = []
		seqs = []
		
		for r in i[1]:
			reads.append(r.split("|")[0])
			seqs.append(r.split("|")[1])
		
		print introns, coverage, ",".join(reads), ",".join(seqs[:4])
		
		
		
		
#		Q = record.letter_annotations["phred_quality"]		
#		upperseq = SeqRecord( Seq(str(record.seq).upper()), id = record.id, description = "" )
#		upperseq.letter_annotations["phred_quality"] = Q
#		print upperseq.format("fastq")
		
			
if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

	
