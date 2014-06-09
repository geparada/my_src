import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq


SeqTable = []

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ..."	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

	f.close()  


def hist_anchor (fastq):
	
	Genome = dict(SeqTable)

	hist_TOTAL = defaultdict(int)
	hist_non_canonical = defaultdict(int)

	print >> sys.stderr, "Contando reads totales desde FASTQ ..."

        f = open(fastq)

        for record in SeqIO.parse(f, "fastq"):
                read = str(record.id)
		introns_sim = read.split('=')[1].split('<>')	
		anchors_totales_FASTQ = read.split('=')[2].split('<>')
		for a in anchors_totales_FASTQ:
			hist_TOTAL[len(record.seq)] += 1


		for i, a in zip(introns_sim, anchors_totales_FASTQ) :

			chr = i.split(':')[0]
			strand = '+'
			if not strand in i:
				strand = '-'
			istart = int(i.split(':')[1].split(strand)[0])
			iend = int(i.split(':')[1].split(strand)[1]) 
			i_dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
			if strand == '-':
				i_dn = i_dn.reverse_complement()
			i_dn = str(i_dn).upper()

		
			if i_dn != 'GTAG' and i_dn != 'GCAG' and i_dn != 'ATAC':
				hist_non_canonical[len(record.seq)] += 1


	for nt in range(1,101):                   #column join
		Total = hist_TOTAL[nt]
		Non_canonical = hist_non_canonical[nt]
	
		print nt, Total, Non_canonical
		
									
			

	f.close()




if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	hist_anchor (sys.argv[2])
