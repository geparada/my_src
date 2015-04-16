import sys
import csv
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
	

def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')

	intron_retention_exon_count = 0
	alt_introns_count = 0
	alt_no_skipper_introns_count = 0
	alt_skipper_introns_count = 0
	alt_exon_variant_introns_count = 0
	shift_count = 0
	non_canonical_shift_count = 0



	for row in reader1:
		
		gene = row[0]
		row = row[1:]
		
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
		intron_retention_exon = row[19].split(",")
		skipped_exons_names = row[20].split(",")
		alt_introns = row[21].split(",")
		alt_no_skipper_introns = row[22].split(",")
		alt_skipper_introns = row[23].split(",")
		alt_exon_variant_introns = row[24].split(",")
		shift = row[25].split(",")
		non_canonical_shift = row[26].split(",")
		
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			
			for i in shift:
				if i!="NO" and (i.split("|")[1]=="GTAG" or i.split("|")[1]=="GCAG" or i.split("|")[1]=="ATAC"):
					pos_donor, pos_aceptor = int(i.split("|")[-2]), int(i.split("|")[-1])
					
					if pos_donor <= 0 and pos_aceptor >= 0:
						micro_exon = Genome[chr][istart+pos_donor:istart] + Genome[chr][iend:iend+pos_aceptor]
						intron_seq = str(Genome[chr][istart:iend]).upper()
						
						if strand == "-":
							micro_exon = str(Genome[chr][iend:iend-pos_donor].reverse_complement() + Genome[chr][istart-pos_aceptor:istart].reverse_complement()).upper()
							intron_seq = str(Genome[chr][istart:iend].reverse_complement()).upper()	  
					
						seq = str("AC" + micro_exon + "AT").upper()

						if seq in intron_seq:
							
							n = 0
							c = 0

							for i1, i2  in zip(intron_seq.split(seq), intron_seq.split(seq)[1:]):
							
								n += len(i1) + 2 + c*len(micro_exon) + c*2
								
								exon_pos = chr + ":" + str(istart+n) + strand + str(istart+n+len(micro_exon))
								
								if strand == "-":
									 exon_pos = chr + ":" + str(iend-n-len(micro_exon)) + strand + str(iend-n)
								
								#print gene, intron, i, i1[-15:]+"AG", micro_exon, "GT"+i2[:15], len(i1), exon_pos
								
								c += 1

					
#					elif pos_aceptor - pos_donor >=3:
#						intron_seq = str(Genome[chr][istart:iend]).upper()
#						up_exon = str(Genome[chr][istart-15:istart]).upper()
#						down_exon = str(Genome[chr][iend:iend+15]).upper()
						
#						if strand == "-":
#							intron_seq = str(Genome[chr][istart:iend].reverse_complement()).upper()	  
#							down_exon = str(Genome[chr][istart-15:istart].reverse_complement()).upper()
#							up_exon = str(Genome[chr][iend:iend+15].reverse_complement()).upper()
						
#						print gene, intron, i, up_exon, intron_seq[:20], intron_seq[-20:], down_exon
							


			
	
					

		
		
			
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2])	
