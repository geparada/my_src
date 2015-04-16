import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from random import randint, sample
from operator import itemgetter
from collections import defaultdict

Transcriptome = {}

Genome = {}
chr_sizes = {}
new_genome = {}
chrs = set([])
splice_sites = defaultdict(list)


def GenomeEditor(chr, micro_exon_istart, micro_exon_iend, strand, new_intron5, new_intron3,  genome=Genome):  #, splice_sites=splice_sites):


	if strand == "+":

		new_chr = Genome[chr][:micro_exon_iend-12] + new_intron3 + Genome[chr][micro_exon_iend:micro_exon_istart] + new_intron5 + Genome[chr][micro_exon_istart+6:]
		Genome[chr] = new_chr		

	if strand == "-":

		new_chr = Genome[chr][:micro_exon_istart-6] + new_intron5 + Genome[chr][micro_exon_istart:micro_exon_iend] + new_intron3 + Genome[chr][micro_exon_iend+12:]
		Genome[chr] = new_chr		

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

	
def Transcriptometabulator(genecode_fasta):
	
	print >> sys.stderr, "Cargando a fasta en la ram ...",
	
	for record in SeqIO.parse(genecode_fasta, "fasta"):
		id = str(record.id).split("|")[0]
		Transcriptome[id] = record.seq
		
	print >> sys.stderr, "OK"
	
		
	
def main(bed12):

	print >> sys.stderr, "Extrayendo intrones del bed12 ...",

         # row[0]  chr
         # row[1]  alignment start
         # row[2]  alignment end
         # row[3]  name 
         # row[4]    
         # row[5]  strand
         # row[6]  aligment start
         # row[7]  aligment end
         # row[8]  
         # row[9] blocknum
         # row[10] blocksizes
         # row[11] qstarts  


	simulation_genome = open("simulation_genome.fa", 'w') 

	csv.field_size_limit(1000000000)

	n = 99
	min_intron_lenght = 80 	

	introns = defaultdict(set)
	acepted_introns = set([])

	transcript_intron_info = {}



	for row in csv.reader(open(bed12), delimiter = '\t'):  #Generando lista de intrones no solapantes
		
		
		qName = row[3]

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]
		qstart = 0

		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
			
			qstart = qstart + b
			tag_start = qstart - n
			tag_end = qstart + n

			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = row[0] + ":" +  str(istart) + row[5] + str(iend)	
			intron = chr + ":" + str(istart) + strand + str(iend)
			ilength = iend - istart

			dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
			intron5 = Genome[chr][istart:istart+6]
			intron3 = Genome[chr][iend-12:iend]
					
			if strand == "-":
				intron3 = Genome[chr][istart:istart+12]
				intron5 = Genome[chr][iend-6:iend]
				dn = dn.reverse_complement()

			dn = str(dn).upper()

			if ilen >= 250 and (dn=="GTAG" or dn=="GCAG" or dn=="ATAC"):
				introns[chr].add((intron, chr, istart, iend))
				splice_sites["intron5" + strand].append(intron5)
				splice_sites["intron3" + strand].append(intron3)



	for c in introns.items():

		chr, intron_list = c

		previous_end = 0

		for i in sorted(intron_list, key=itemgetter(2)):

			intron, chr, istart, iend = i

			if istart > previous_end:

				acepted_introns.add(intron)
				previous_end = iend



	for row in csv.reader(open(bed12), delimiter = '\t'):
		
		try:
		
			qName = row[3]
			seq = Transcriptome[qName]

			qstarts = map (int, row[11].strip(",").split(","))                      
			blocksizes = map(int, row[10].strip(",").split(","))

			start = int(row[1])
			strand = row[5]
			bn = int(row[9])
			chr = row[0]
			qstart = 0

			for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
				
				qstart = qstart + b
				tag_start = qstart - n
				tag_end = qstart + n

				istart = start + q1 + b
				iend = start + q2
				ilen = iend - istart
				intron = row[0] + ":" +  str(istart) + row[5] + str(iend)	
				intron = chr + ":" + str(istart) + strand + str(iend)
				ilength = iend - istart
				intron_seq = Genome[chr][istart:iend]
				
				if strand == '+' :                          #Para los que aliniean en la hebra +
								   
					if tag_start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
						tag_start = 0  
					if tag_end>len(seq):
						tag_end=len(seq)
					tag = seq[tag_start:tag_end]
								  
				if strand == '-' :

					intron_seq = intron_seq.reverse_complement()
				
					if tag_end>len(seq):                 #Para los que alinian en la hebra - es todo al inverso
						tag_end=len(seq)                                       
					tag = seq[-tag_end:-tag_start]
					if tag_start<=0: 
						tag = seq[-tag_end:]
										 
				if len(tag) == 2*n: #and (intron in acepted_introns):

					info = intron_seq, min_intron_lenght, tag, chr, istart, iend, strand
					transcript_intron_info[intron] = info


		except KeyError:
			pass


	for intron in acepted_introns:


		try:

			intron_seq, min_intron_lenght, tag, chr, istart, iend, strand = transcript_intron_info[intron]


			core_intron = intron_seq[min_intron_lenght:-min_intron_lenght]
			micro_exon_start = randint(0, len(core_intron))
			micro_exon_lenght = randint(1,25)
			micro_exon_tag = tag[:n] + core_intron[micro_exon_start:micro_exon_start + micro_exon_lenght] + tag[n:]

			micro_exon_iend = istart + min_intron_lenght + micro_exon_start
			micro_exon_istart = micro_exon_iend + micro_exon_lenght

			if strand == "-":
				micro_exon_istart = iend - (min_intron_lenght + micro_exon_start + micro_exon_lenght)  #Los di vuelta
				micro_exon_iend = micro_exon_istart + micro_exon_lenght


			intron5 = Genome[chr][micro_exon_istart:micro_exon_istart+6]
			intron3 = Genome[chr][micro_exon_iend-12:micro_exon_iend]
						
			if strand == "-":
				intron3 = Genome[chr][micro_exon_iend:micro_exon_iend+12].reverse_complement()         #cambiado
				intron5 = Genome[chr][micro_exon_istart-6:micro_exon_istart].reverse_complement()

			new_intron5 = sample(splice_sites["intron5" + strand],1)[0]
			new_intron3 = sample(splice_sites["intron3" + strand],1)[0]

			GenomeEditor(chr, micro_exon_istart, micro_exon_iend, strand, new_intron5, new_intron3)


			print intron, micro_exon_tag, core_intron[micro_exon_start:micro_exon_start + micro_exon_lenght], micro_exon_start, micro_exon_start + micro_exon_lenght, intron3, intron5, micro_exon_iend, micro_exon_istart, new_intron3, new_intron5

		except KeyError:
			pass


	for i in Genome.items():

		chr, seq = i
		seq = str(seq)
		simulation_genome.write(">" + chr + "\n")
		simulation_genome.write(seq + "\n") 	


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	Transcriptometabulator(sys.argv[2])
	main (sys.argv[3])


#python my_src/Micro-exons/micro_exon_simulation.py db/genomes/hg19.fa db/trasncriptomes/Gene_models/gencode/v19/gencode.v19.pc_transcripts.fa db/trasncriptomes/Gene_models/gencode/v19/gencode.v19.annotation.bed

#						if strand == "+":
#
#							for n, c in zip(new_intron5, range(micro_exon_istart, micro_exon_istart+6)):
#								new_genome[chr + ":" + c] = n
#
#							for n, c in zip(new_intron3, range(micro_exon_iend-12, micro_exon_iend)):
#								new_genome[chr + ":" + c] = n
#
#						if strand == "-":
#
#							for n, c in zip(new_intron3, range(micro_exon_istart, micro_exon_istart+12)):
#								new_genome[chr + ":" + c] = n
#
#							for n, c in zip(new_intron5, range(micro_exon_iend-6, micro_exon_iend)):
#								new_genome[chr + ":" + c] = n	