import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

hexamers = set([])
ESE53 = {}
ESE = {}
ESS = {}
ISE = {}


def hexamer_reader (hexamer_file, hexamer_dict):
	""" reads the haxamers file """
	
	hexamer_groups = defaultdict(list)
	
	for row in csv.reader(open(hexamer_file), delimiter = '\t'):
		seq, group = row
		hexamers.add(seq)
		hexamer_groups[group].append(seq)
	
	for i in hexamer_groups.items():
		group, seqs = i
		
		c = 0		
		for seq in seqs:
			c += 1
			ID = group + "_" + str(c)
			hexamer_dict[seq] = ID
				

def main (ESE5_file, ESE3_file, ESS_file, ISE_file, hg19):
	
	""" Generates a BED file with the hits positions """

	hexamer_reader(ESE3_file, ESE)	
	hexamer_reader(ESE5_file, ESE)
	hexamer_reader(ESS_file, ESS)	
	hexamer_reader(ISE_file, ISE)
	
	ESE_uniq = defaultdict(list)
	
	Genome = []
		
	f = open(hg19)
	for chrfa in SeqIO.parse(f, "fasta"):
		chr = chrfa.id
		seq = chrfa.seq
		Genome.append((chr,seq))
			
	for chrfa in Genome:
		chr, seq = chrfa
		kmer_len = 6
		strand = "+"
		blockCount = "1"
		blockSizes = str(kmer_len)
		blockStarts = "0"
		
		for i in range(len(seq)- kmer_len):
			start = i
			end = i + kmer_len 
			chr_hex = str(seq[start:end]).upper()
			
			try:
				print "\t".join([chr, str(start), str(end ), ESE[chr_hex], "0", strand, str(start), str(end), "200,0,0", blockCount, blockSizes, blockStarts])
				
			except KeyError:
				pass

			try:		
				print "\t".join([chr, str(start), str(end ), ESS[chr_hex], "0", strand, str(start), str(end), "0,200,0", blockCount, blockSizes, blockStarts])
				
			except KeyError:
				pass
							
			try:		
				print "\t".join([chr, str(start), str(end ), ISE[chr_hex], "0", strand, str(start), str(end), "0,0,200", blockCount, blockSizes, blockStarts])
	
			except KeyError:
				pass
		
		strand = "-"
		seq = seq.reverse_complement()
				
		for i in range(len(seq)- kmer_len):
			start = i
			end = i + kmer_len 
			chr_hex = str(seq[start:end]).upper()
			
			try:
				print "\t".join([chr, str(start), str(end ), ESE[chr_hex], "0", strand, str(start), str(end), "200,0,0", blockCount, blockSizes, blockStarts])
				
			except KeyError:
				pass

			try:		
				print "\t".join([chr, str(start), str(end ), ESS[chr_hex], "0", strand, str(start), str(end), "0,200,0", blockCount, blockSizes, blockStarts])
				
			except KeyError:
				pass
							
			try:		
				print "\t".join([chr, str(start), str(end ), ISE[chr_hex], "0", strand, str(start), str(end), "0,0,200", blockCount, blockSizes, blockStarts])
	
			except KeyError:
				pass	
	
	
			



if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])	
