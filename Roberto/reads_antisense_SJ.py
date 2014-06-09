import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
from Bio.SeqRecord import SeqRecord


read_introns = defaultdict(set)
antisense = defaultdict(list)

def main(bam):

	samfile = pysam.Samfile(bam)
	

	
	for read in samfile.fetch():
		
		if read.flag!=4:
			
			qname = read.qname
			intron = samfile.getrname(read.tid).split("|")[1]
			flag = read.flag
			tstart =  read.pos
			matchs = read.cigar[0][1]
			tag_len = 100
			seq = read.seq
			qual = read.qual
			
			if (tag_len/2 - 8) >= tstart and (tag_len/2 + 8) <= tstart + matchs:
			
				read_introns[qname].add((intron, flag))

				if flag == 16:
					antisense[qname].append((intron, flag, tstart, matchs, seq, qual))
	
	sense_count = 0
	antisense_count = 0 
	
	for i in read_introns.items():
		introns = list(i[1])[0]
		if len(read_introns[qname])==1:
			flag = introns[1]
			if flag == 0:
				sense_count += 1
			elif flag == 16:
				antisense_count += 1
	
#	print sense_count, antisense_count

	intron_reads_antisense = defaultdict(list)
			
	for i in antisense.items():
		qname = i[0]
		intron = i[1][0][0]
		flag = i[1][0][1]
		tstart = i[1][0][2]
		matchs = i[1][0][3]
		seq = i[1][0][4]
		qual = i[1][0][5]
#		print qname, intron, flag, tstart, matchs, (read_introns[qname])
		
		if len(read_introns[qname])==1:
			print "@" + qname
			print seq
			print "+"
			print qual

#			fastq = SeqRecord( seq, id = qname, description = "" )
#			fastq.letter_annotations["phred_quality"] = qual
		
		
#			print fastq.format("fastq"),	

	
#			intron_reads_antisense[intron].append((qname, seq))

			

	
#	for i in intron_reads_antisense.items():
#		qnames = []
#		seqs = []
#		for n in range(len(i[1])):
#			qnames.append(i[1][n][0])
#			seqs.append(i[1][n][1])
		
		#print i[0], len(i[1]), ",".join(qnames), ",".join(seqs[:5]) 
	
			
			
		

#Si es que alinea a dos splice juntions, no sirve				
			
			
	
		


if __name__ == '__main__':
	main(sys.argv[1])
