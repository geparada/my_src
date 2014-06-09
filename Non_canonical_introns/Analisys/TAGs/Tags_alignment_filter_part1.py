import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


TAG_seq = {}
fasta = open("reads_tags.filter.fa", "w")

	
def TAGstabulator(tags_fasta):
	f = open(tags_fasta)
	for chrfa in SeqIO.parse(f, "fasta"):
		TAG_seq[str(chrfa.id)] = str(chrfa.seq)	
	f.close()



def main(reads_genome_bam, reads_tags_SAM):
	
	reader1 = csv.reader(open(reads_tags_SAM), delimiter = '\t')

	samfile1 = pysam.Samfile( reads_genome_bam, "rb" )
	samfile2 = pysam.Samfile( reads_genome_bam, "rb" )
	
	spliced_and_multimapping = {}
	
		
	for row in samfile1:
		
		
		try:
			read = row.qname
			multimapping_tag = row.tags[2][1]
			for segment in row.cigar:               #Cuando se encuentra con un 3, significa que hay un intron en la linea CIGAR
				if segment[0] == 3:
					spliced_and_multimapping[read] = row.cigar 
			if multimapping_tag !=1:                #Revisa los tags generados por MapSplice y verifica si hace multimapping
				spliced_and_multimapping[read] = row.tags[2]
			
		except IndexError:
			pass
	
	introns_finded = {}
	

	
	for row in reader1:
		read = row[0]
		flag = row[1]
		tag_name = row[2]
		start = int(row[3])           #Sam es 1 referenciado 
		cigar = row[5]
		seq = row[9]
		qual = row[10]
		intron_tag = tag_name.split("|")[0]		
		
		mismatches = int(row[13].strip("NM:i"))
		matches = int(cigar.strip('M'))
		fist_block_tag = int(tag_name.split("|")[2].split("_")[0])
		anchor_up =  fist_block_tag - start
		anchor_down =  start + matches - fist_block_tag

		
		if mismatches <= 2 and anchor_up >= 8 and anchor_down >= 8 and spliced_and_multimapping.has_key(read)==False:
			read_window = seq[fist_block_tag-start-8+1:fist_block_tag-start+8+1]
			seq_tag = TAG_seq[tag_name]
			tag_window = seq_tag[fist_block_tag-8:fist_block_tag+8]
			if read_window == tag_window:
				introns_finded[read] = (seq, qual, intron_tag, anchor_up, anchor_down)
	
#	for i in introns_finded:
#		read = i[0]
#		seq = i[1]
#		qual = i[2]
#		intron_tag = i[3]
#		anchor_up = i[4]
#		anchor_down = i[5]
		
	for row in samfile2:
		read = row.qname
		try:
			seq = introns_finded[read][0]
			qual = introns_finded[read][1]
			intron_tag = introns_finded[read][2]
			anchor_up = introns_finded[read][3]
			anchor_down = introns_finded[read][4]
			
			print read, seq, qual, intron_tag, anchor_up, anchor_down, row.is_unmapped
			
			if row.is_unmapped:
				print read, seq, qual, intron_tag, anchor_up, anchor_down
				fasta.write(">" + read + "\n")
				fasta.write(seq + "\n")				
							
			elif (row.qlen - (row.qend - row.qstart) ) > 8:            #Si el read es mapeado... solo se admite si deja 8 nucleotidos sin mapear en algun extremo
				print read, seq, qual, intron_tag, anchor_up, anchor_down
				fasta.write(">" + read + "\n")
				fasta.write(seq + "\n")
									
		except KeyError:
			pass


#Este script filtra todos aquellos tags que fueron mapeados anteriormente por MapSplice a un sitio de splicing o a una zona del genoma.
#Solo se admitio el read cuando hubo al menos 8 nucleotidos que hicieron softcpling		


if __name__ == '__main__':
	TAGstabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3])
