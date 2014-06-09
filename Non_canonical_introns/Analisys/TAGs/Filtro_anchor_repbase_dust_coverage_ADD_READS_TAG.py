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

def main(SJ_introns, repeats, dust_IDs, min_coverage):
	Genome = dict(SeqTable)
	repeat_list = []
	read_seq_list = []
	intron_read_seq = defaultdict(set)
	intron_info = defaultdict(set)
	
	
	repeat_list = []
	read_seq_list = []
	
	intron_read_seq = defaultdict(set)
	intron_info = defaultdict(set)
	
	reader1 = csv.reader(open(SJ_introns), delimiter = ' ')
	reader2 = csv.reader(open(repeats), delimiter =  '\t')
	reader3 = csv.reader(open(dust_IDs), delimiter =  ' ')
				
	
	print >> sys.stderr, "Cargando reads con repeticiones en la ram ...",
	for row in reader2:
				
		matches = int(row[0])
		qName = row[9]
		qSize = int(row[10])
		rep_percent = percent(matches, qSize)
		
		if rep_percent > 50:
			repeat_list.append((qName, rep_percent))
			
	for row in reader3:
		ID = row[0][1:]
		repeat_list.append((ID, 'DUST'))
		
	repeat_dict = dict(repeat_list)
	print >> sys.stderr, "OK"

	print >> sys.stderr, "Aplicando filtros...",
	for row in reader1:

		read = row[0]
		chr = row[1]
		istart = int(row[2])
		iend = int(row[3])
		strand = row[4]
		ilen = int(row[5])
		intron = row[6]
		dn = row[7]
		start = int(row[8])
		cigar = row[9]
		e5s = int(row[10])
		e5e = int(row[11])
		e3s = int(row[12])
		e3e = int(row[13])
		seq = row[14]
		end = e3e
		
		if dn=='GTAG' or dn=='GCAG' or dn=='ATAC':
	
			intron_read_seq[intron].add(seq)                                #diccionario intron - secuecias distintas que lo apollan			
			intron_info[intron].add((chr, istart, iend, strand, ilen))      #diccionario intro - info
		
		elif repeat_dict.has_key(read)==False:                            #El filtro de repeticiones solo se aplica a los no canonicos

			intron_read_seq[intron].add(seq)                                #diccionario intron - secuecias distintas que lo apollan			
			intron_info[intron].add((chr, istart, iend, strand, ilen))      #diccionario intro - info			
			
					
	
	#intron_info_dict = dict(intron_info_list)				
			 
	for pairs in intron_read_seq.items():
		intron = pairs[0]
		coverage = len(pairs[1])
		
		if coverage >= min_coverage:          #Aplicando filtro de coverage
			Introns.append((intron, coverage, list(intron_info[intron])))
	print >> sys.stderr, "OK"



def no_can_reads(SJ_introns, Genome):
	non_can_reads_list = defaultdict(list)
	reader1 = csv.reader(open(SJ_introns), delimiter = ' ')
	
	for row in reader1:

		read = row[0]
		chr = row[1]
		istart = int(row[2])
		iend = int(row[3])
		strand = row[4]
		ilen = int(row[5])
		intron = row[6]
		dn = row[7]
		start = int(row[8])
		cigar = row[9]
		e5s = int(row[10])
		e5e = int(row[11])
		e3s = int(row[12])
		e3e = int(row[13])
		seq = row[14]
		end = e3e
				
		if dn!='GTAG' and dn!='GCAG' and dn!='ATAC':
			non_can_reads_list[intron].append(read)
	
	return non_can_reads_list



def dinucleotides(SJ_introns):
	Genome = dict(SeqTable)
	
	print >> sys.stderr, "Extrayendo reads de intrones no canonicos ...",
	
	non_can_reads_list = no_can_reads(SJ_introns, Genome)
	non_can_reads_dict = dict(non_can_reads_list)
	
	print >> sys.stderr, "OK"
	
	for row in Introns:               #La estructura es del el tipo ('chr2:102972860+102979092', 11, [('chr2', 102972860, 102979092, "+", 6232)])
		
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

		if non_can_reads_dict.has_key(intron):
			
			reads = ",".join(non_can_reads_dict[intron])
						
			print chr, istart, iend,  intron + '|' + str(dn).upper() + '|' + reads, coverage, strand         #Este formato es necesario para hacer un liftover

		else:

			print chr, istart, iend,  intron + '|' + str(dn).upper() + '|' + "0", coverage, strand
			
			#Revisar la estrutura de la lista donde se esta iterando, ya que pasa algo raro con la variable ilen ...
			
if __name__ == '__main__':
	
	main(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])) 
	Genomictabulator(sys.argv[1])
	dinucleotides(sys.argv[2])
	
	#USAGE: python Filtro_anchor_repbase_coverage.py genome.fa SJ_introns repeats_repbase.psl low_complexity_IDs.dust min_coverage  
