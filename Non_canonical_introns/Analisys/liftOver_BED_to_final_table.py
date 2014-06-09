import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable = []
		
def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

	f.close() 

def main(BED):
	Genome = dict(SeqTable)
	reader = csv.reader(open(BED), delimiter = '\t')
	
	for row in reader:
		chr = row[0]
		istart = int(row[1])
		iend = int(row[2])
		intron_mm9 = row[3].split("|")[0]
		dn_mm9 = row[3].split("|")[1]
		cDNA_coverage = int(row[3].split("|")[2])
		EST_coverage = int(row[3].split("|")[3])
		cDNA = row[3].split("|")[4]
		EST = row[3].split("|")[5]
		strand = row[5]
		
		intron = chr + ':' + str(istart) + strand + str(iend)
		ilenght = iend - istart
		
		dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]

		if strand == '-':
			dn = dn.reverse_complement()
			
		dn = str(dn).upper()
			
		if dn_mm9 == dn:
			if dn != "GTAG" and dn != "GCAG" and dn != "ATAC":
				print intron, cDNA_coverage + EST_coverage, chr, strand, istart, iend, ilenght, dn, cDNA + ',' + EST
			
			else:
				print intron, cDNA_coverage + EST_coverage, chr, strand, istart, iend, ilenght, dn, 0				








if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])
