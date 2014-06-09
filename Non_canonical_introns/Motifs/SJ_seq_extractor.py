import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()


def main (final_table):
	""" Extrae secuencias de SJ para analisis de motivos """
	
	reader = csv.reader(open(final_table), delimiter = '\t')
	
	can_1 = open("canonicos_exon5", 'w') 
	can_2 = open("canonicos_intron5", 'w')
	can_3 = open("canonicos_intron3", 'w') 
	can_4 = open("canonicos_exon3", 'w')
	
	no_can5_1 = open("no_canonicos5_exon5", "w")
	no_can5_2 = open("no_canonicos5_intron5", "w")
	no_can5_3 = open("no_canonicos5_intron3", "w")	
	no_can5_4 = open("no_canonicos5_exon3", "w")
	
	no_can3_1 = open("no_canonicos3_exon5", "w")
	no_can3_2 = open("no_canonicos3_intron5", "w")
	no_can3_3 = open("no_canonicos3_intron3", "w")	
	no_can3_4 = open("no_canonicos3_exon3", "w")	
	
	for row in reader:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = row[8]
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		bodymap_seq = row[16]
		gm12878_seq = row[17]
		DR = row[18]
		
		L = 300
				
		ext_seq = L
		in_seq = L
		
		if in_seq > ilength:
			in_seq = ilength
		
		if dn=='GTAG' or dn=='GCAG' or dn=='ATAC':
			exon5 = Genome[chr][istart-300:istart]
			intron5 = Genome[chr][istart:istart+300]
			intron3 = Genome[chr][iend-300:iend]
			exon3 = Genome[chr][iend:iend+300]		
			
			if strand == "-":
				exon3 = Genome[chr][istart-300:istart].reverse_complement()
				intron3 = Genome[chr][istart:istart+300].reverse_complement()
				intron5 = Genome[chr][iend-300:iend].reverse_complement()
				exon5 = Genome[chr][iend:iend+300].reverse_complement()				 
			
			exon5 = str(exon5).upper() + "\n"
			intron5 = str(intron5).upper() + "\n"
			intron3 = str(intron3).upper() + "\n"
			exon3 = str(exon3).upper() + "\n"
			
			
			can_1.write(exon5)
			can_2.write(intron5)
			can_3.write(intron3)
			can_4.write(exon3)
		
		elif dn[2:] == "AG" or dn[2:] == "AC"
			
			
				
	can_1.close()
	can_2.close()	
	can_3.close()
	can_4.close()

#Desafio: Las secuencias exonicas deben ser exclusibamente exonicas y las intrincas, exclucibamente intronicas

	
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2]) 


