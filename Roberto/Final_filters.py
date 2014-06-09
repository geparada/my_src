import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict

Genome = {}


def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def Genomictabulator(fasta):

	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def main(bam):

	samfile = pysam.Samfile(bam, "rb")
	out = pysam.Samfile(sys.argv[2].strip(".bam") + ".final_filters.bam", "wb", template=samfile)
	SJ_NM_0 = defaultdict(int)
	SJ_NM_x = defaultdict(int)
	black_list = set([]) 

	for read in samfile.fetch():
		qname, start, end, seq_len, mismatches, flag, chr, gene, intron = read.qname, read.pos, read.aend, read.qlen, read.tags[4][1], read.flag, samfile.getrname(read.tid), read.tags[7][1], read.tags[8][1]
		
		if seq_len < 25:
			if mismatches==0:
				SJ_NM_0[intron] += 1
			else:
				SJ_NM_x[intron] += 1	
	
	for i in SJ_NM_x.items():
		intron, NM_x_count = i
		NM_0_count = SJ_NM_0[intron]
				
		if percent(NM_x_count, NM_x_count + NM_0_count) >= 75:
			black_list.add(intron)

		
	for read in samfile.fetch():
		qname, start, end, seq_len, mismatches, flag, chr, gene, intron = read.qname, read.pos, read.aend, read.qlen, read.tags[4][1], read.flag, samfile.getrname(read.tid), read.tags[7][1], read.tags[8][1]
		strand = "+"
		if read.is_reverse:
			strand = "-"

		NN = str(Genome[chr][start-2:start])
		if strand=="-":
			NN = str(Genome[chr][end:end+2].reverse_complement())		

		if NN!="CG" and (intron in black_list)==False:
			out.write(read)





			
	samfile.close()

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2])

