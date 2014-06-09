import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[str(chrfa.id)] = chrfa.seq

	f.close()

def main(reads_tags):
	reader1 = csv.reader(open(reads_tags), delimiter = ' ')	
	
	for row in reader1:
		read = row[0]
		seq = row[1]
		qual = row[2]
		intron = row[3]
		up_anchor = int(row[4])
		down_anchor = int(row[5])
		
		chr = intron.split(":")[0]
		strand = ""
		if "-" in intron:
			strand = "-"
		elif "+" in intron:
			strand = "+"
		
		istart = int(intron.split(":")[1].split(strand)[0]) 
		iend = int(intron.split(":")[1].split(strand)[1])
		
		ilen = iend - istart



		dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
		if strand == '-':
			dn = dn.reverse_complement()
			
		dn = str(dn).upper()

		start = int(0)
		cigar = "FINDED_BY_TAG"
		e5s = int(0)
		e5e = int(0)
		e3s = int(0)
		e3e = int()
		end = e3e
		
		print read, chr, istart, iend, strand, ilen, intron, dn, start, cigar, e5s, e5e, e3s, e3e, seq, end  


	

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])
