import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC


Genome = {}


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[str(chrfa.id)] = chrfa.seq

	f.close()  

	print >> sys.stderr, "OK"

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')
	
	derivados_GC = defaultdict(int)
	indef_GC = defaultdict(int)
	canonicos_GC = defaultdict(int)
	
	derivados = 0
	indef = 0
	canonicos = 0			

	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)

#	Canonicos = []
	
#	U2_GTAG = []
#	U12_ATAC = []
#	U12_GTAG = []
#	indef = []

	for row in reader1:

		gene = row[0]
		row = row[1:]
		
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = float(row[8])
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		bodymap_seq = row[16]
		gm12878_seq = row[17]
		DR = row[18]
		intron_retention_exon = row[19].split(",")
		skipped_exons_names = row[20].split(",")
		alt_introns = row[21].split(",")
		alt_no_skipper_introns = row[22].split(",")
		alt_skipper_introns = row[23].split(",")
		alt_exon_variant_introns = row[24].split(",")
		shift = row[25].split(",")
		non_canonical_shift = row[26].split(",")

		intron_seq = Genome[chr][istart:iend]
		if strand == "-":
			intron_seq = intron_seq.reverse_complement()
		GC_int = int(str(GC(intron_seq)/100)[2]) + 10*int(str(GC(intron_seq)/100)[0])
		
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				derivados_GC[GC_int] += 1
				derivados += 1		
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				derivados_GC[GC_int] += 1		
				derivados += 1
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				derivados_GC[GC_int] += 1
				derivados += 1				
					
			elif shift[0]=="NO":
				indef_GC[GC_int] += 1
				indef += 1
		
		else:
			
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG :
				canonicos_GC[GC_int] += 1
				canonicos += 1
					
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC:
				canonicos_GC[GC_int] += 1
				canonicos += 1		
					
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG:
				canonicos_GC[GC_int] += 1				
				canonicos += 1		
	
	print "indefinidos", "derivados", "canonicos"
		
	for i in range(11):
		print i*10, percent(indef_GC[i], indef) , percent(derivados_GC[i], derivados), percent(canonicos_GC[i], canonicos) 

#	print 100, percent(indef_GC[10], indef) , percent(derivados_GC[10], derivados), percent(canonicos_GC[i], canonicos) 			
					

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])

#	Genomictabulator("/home/geparada/db/genome/hg19.fa")
#	main("../TOTAL/FINAL_TABLE.alternative_splicing.fold.highfold.gene_names")	
	
	
#~/db/genome/hg19.fa ../TOTAL/FINAL_TABLE.alternative_splicing.fold.highfold.gene_names
