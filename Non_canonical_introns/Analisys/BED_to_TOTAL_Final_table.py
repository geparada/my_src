import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def main(final_table):
	reader1 = csv.reader(open(final_table), delimiter = '\t')
	
	for row in reader1:
		
		chr = row[0]
		istart = int(row[1])
		iend = int(row[2])
		strand = row[5]
		
		dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
			
		if strand == '-':
			dn = dn.reverse_complement()
				
		dn = str(dn).upper()
		
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			print row[3].split("|")[0], chr, strand, istart, iend, " ".join(row[3].split("|")[5:])
		
		#print row[3].replace("|", " ")
		

		




if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2]) 
