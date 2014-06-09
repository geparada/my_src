import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}
		
def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	f.close() 
	

def main(bed):
	csv.field_size_limit(1000000000)
	reader = csv.reader(open(bed), delimiter = '\t')
	
	for row in reader:
		chr = row[0]
		istart = int(row[1])
		iend = int(row[2])
		intron = row[3].split("|")[0]
		dn_type = row[3].split("|")[1]
		dn_type_score = row[3].split("|")[2]
		reads = row[3].split("|")[3]
		coverage = row[4]
		strand = row[5]
		dn = row[8]
		
		human_dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]

		if strand == '-':
			human_dn = human_dn.reverse_complement()
		
		human_dn = str(human_dn).upper()
		human_intron = chr +  ":" + str(istart) + strand + str(iend)
		
		ilength = iend - istart
		
		if human_dn == dn:
			print human_intron, coverage, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, reads
		





if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])
