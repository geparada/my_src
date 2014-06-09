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

SeqTable = []
Introns = []

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando el genoma en la RAM ...",
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = (str(chrfa.id), chrfa.seq)
		SeqTable.append(table)

	f.close()
	print >> sys.stderr, "OK"

def main(bed12_seq, repeats, min_anchor, min_coverage):
	Genome = dict(SeqTable)
	
	repeat_list = []
	read_seq_list = []
	
	intron_read_seq = defaultdict(set)
	intron_info = defaultdict(set)
	
	reader1 = csv.reader(open(bed12_seq), delimiter = '\t')
	reader2 = csv.reader(open(repeats), delimiter =  '\t')
				
	
	print >> sys.stderr, "Cargando reads con repeticiones en la ram ...",
	for row in reader2:
				
		matches = int(row[0])
		qName = row[9]
		qSize = int(row[10])
		rep_percent = percent(matches, qSize)
		
		if rep_percent > 50:
			repeat_list.append((qName, rep_percent))
	
	
	repeat_dict = dict(repeat_list)
	print >> sys.stderr, "OK"

	print >> sys.stderr, "Aplicando filtros...",
	for row in reader1:
		
		chr = row[0]
		alignment_start = row[1]  
		alignment_end = row[2]  
		read = row[3]   
		strand = row[5]  
		start = int(row[6])  
		end = int(row[7])  
		blocknum = map(int, row[9].strip(',').split(','))
		blocksizes = map(int, row[10].strip(',').split(','))
		qstarts = map(int, row[11].strip(',').split(','))
		seq = row[12]
		
		if repeat_dict.has_key(read)==False:                          #Filtrando los repeats
		
			for q1, q2, b1, b2 in zip(qstarts, qstarts[1:], blocksizes, blocksizes[1:]):
				istart = start + q1 + b1
				iend = start + q2
				ilen = iend - istart
				intron = row[0] + ':' +  str(istart) + row[5] + str(iend)
										
				if b1 >= min_anchor and b2 >= min_anchor and ilen >= 40:          #Filtrando anchor y length

						
					intron_read_seq[intron].add(seq)    #diccionario intron - secuecias distintas que lo apollan			
					#intron_reads[intron].append(read)                   #diccionario intron - reads
					intron_info[intron].add((chr, istart, iend, strand, ilen))
	
	#intron_info_dict = dict(intron_info_list)				
			 
	for pairs in intron_read_seq.items():
		intron = pairs[0]
		coverage = len(pairs[1])
		
		if coverage >= min_coverage:
			Introns.append((intron, coverage, list(intron_info[intron])))
	print >> sys.stderr, "OK"

	repeat_list = []
	read_seq_list = []
	intron_read_seq = defaultdict(set)
	intron_info = defaultdict(set)


def dinucleotides():
	Genome = dict(SeqTable)
	for row in Introns:               #La estructura es del el tipo ('chr2:102972860+102979092', 11, [('chr2', 102972860, 102979092, 6232)])
		#print row
		intron = row[0]
		coverage = row[1]
		chr = row[2][0][0]                   
		istart = row[2][0][1]
		iend = row[2][0][2]
		strand = row[2][0][3]
		ilen = row[2][0][4]
		
		dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
		if strand == '-':
			dn = dn.reverse_complement()
			
		print intron, coverage, chr, strand, istart, iend, ilen, str(dn).upper()



if __name__ == '__main__':
	
	main(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])) 
	Genomictabulator(sys.argv[1])
	dinucleotides()
	
	#USAGE: python Filtro_anchor_repbase_coverage.py genome.fa bed12_seq repeats min_anchor min_coverage  
