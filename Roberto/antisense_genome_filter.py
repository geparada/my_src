import sys
import csv
import pysam
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

refseq_gene_ID = {}
read_introns = defaultdict(set)
antisense = defaultdict(list)

def Get_gene_ID(ID_list):
	
	for row in csv.reader(open(ID_list), delimiter = '\t' ):
		refseq_gene_ID[row[0]] = row[1]
		

def main(sam_tags, sam_genome):
	
	tags = pysam.Samfile(sam_tags)
	genome = pysam.Samfile(sam_genome)
	
	reads_genome_mismatch = {}
		
	for read in genome.fetch():
		
		if read.flag!=4:
			
			qname = read.qname
			chr = genome.getrname(read.tid)
			flag = read.flag
			tstart =  read.pos
			matchs = read.cigar[0][1]
			tag_len = 100
			seq = read.seq
			qual = read.qual
			mismatch = read.tags[2][1]
			
			reads_genome_mismatch[qname] = mismatch
	
	for read in tags.fetch():
		if read.flag!=4:
			
			qname = read.qname
			refseq = tags.getrname(read.tid).split("|")[0]
			intron = tags.getrname(read.tid).split("|")[1]
			flag = read.flag
			tstart =  read.pos
			matchs = read.cigar[0][1]
			tag_len = 100
			seq = read.seq
			qual = read.qual
			mismatch = read.tags[2][1]

			
			if (tag_len/2 - 8) >= tstart and (tag_len/2 + 8) <= tstart + matchs:
			
				read_introns[qname].add((intron, flag))

				if flag == 16:
					antisense[qname].append((refseq, intron, flag, tstart, matchs, mismatch, seq, qual))
					

	for i in antisense.items():
		qname = i[0]
		refseq = i[1][0][0]
		intron = i[1][0][1]
		flag = i[1][0][2]
		tstart = i[1][0][3]
		matchs = i[1][0][4]
		mismatch_tags = i[1][0][5]
		seq = i[1][0][6]
		qual = i[1][0][7]
		gene = refseq_gene_ID[refseq]
		
		mismatch_genome = 0
		if reads_genome_mismatch.has_key(qname):
			mismatch_genome = reads_genome_mismatch[qname]	
		
		if len(read_introns[qname])==1:
			if mismatch_tags < mismatch_genome:
				print qname, seq, gene, refseq, intron, "tags_better_genome", mismatch_tags, mismatch_genome 				
				
			if mismatch_tags == mismatch_genome:			
				print qname, seq, gene, refseq, intron, "tags_equal_genome", mismatch_tags, mismatch_genome 	

			if mismatch_tags > mismatch_genome:			
				print qname, seq, gene, refseq, intron, "genome_better_tags", mismatch_tags, mismatch_genome 	


#Falta  ponerle nombre de los genes


if __name__ == '__main__':
	Get_gene_ID(sys.argv[1])
	main(sys.argv[2], sys.argv[3])
