import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict


Genome = {}
blocks_starts = set([])
blocks_ends = set([])

donor_aceptor_coverge = defaultdict(int)
exons_junctions = defaultdict(set)

score_U2_GTAG = float(63.2891027914)
score_U12_ATAC = float(60.9280810964)
score_U12_GTAG = float(61.4553595446)

	

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()
	
	
def Exonextractor (gencode_BED12, final_table):
	""" Genera dicionarios de coodenadas de exones que estan presentes en gencode """

	for row in csv.reader(open(final_table), delimiter = '\t'):
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = int(row[5])
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
		tissues_coverage = int(row[16])
		tissues_number = row[17]
		tissues_names = row[18]

		
		total_coverage = bodymap_coverage + gm12878_coverage + tissues_coverage
		
		donor_aceptor_coverge[chr + strand + str(istart)] += total_coverage
		donor_aceptor_coverge[chr + strand + str(iend)] += total_coverage			
	
	
	for row in csv.reader(open(gencode_BED12), delimiter = '\t'):
		chr = row[0]
		start = int(row[1])
		end = row[2]
		name = row[3]
		strand = row[5]
		block_sizes = map (int, row[10].strip(",").split(","))		
		tstarts_0 = map (int, row[11].strip(",").split(","))
		
		
		for t, b in zip(tstarts_0[1:-1], block_sizes[1:-1]):      #No se toman encuenta exones terminales ya que hay muchas isoformas con 5' o 3' alternativos que generan interrupciones de los tags
			
			block_start = chr + strand + str(start + t)
			block_end = chr + strand + str(start + t + b)
			
			blocks_starts.add(block_start)
			blocks_ends.add(block_end)
			
			exons_junctions[block_start].add(((chr, strand, str(start + t + b)), donor_aceptor_coverge.get(block_end, 0)))
			exons_junctions[block_end].add(((chr, strand, str(start + t))  , donor_aceptor_coverge.get(block_start, 0)))

	

def main (final_table, L):
	""" Extrae secuencias de SJ para analisis de motivos
	La idea no es mezclar secuenencias exonicas con intronicas
	Los tags se interrumpen si es que se encuentran con otro splice juntion, de manera que siempre se escoge el exon mas corto"""
	
	
	can_1 = open("canonicos_exon5", 'w') 
	can_2 = open("canonicos_intron5", 'w')
	can_3 = open("canonicos_intron3", 'w') 
	can_4 = open("canonicos_exon3", 'w')
	can_5 = open("canonicos_intron_up", 'w') 
	can_6 = open("canonicos_exon_up", 'w')
	can_7 = open("canonicos_exon_down_", 'w') 
	can_8 = open("canonicos_intron_down", 'w')	
	
	alt_can_1 = open("canonicos_alternativos_exon5", 'w') 
	alt_can_2 = open("canonicos_alternativos_intron5", 'w')
	alt_can_3 = open("canonicos_alternativos_intron3", 'w') 
	alt_can_4 = open("canonicos_alternativos_exon3", 'w')
	alt_can_5 = open("canonicos_alternativos__intron_up", 'w') 
	alt_can_6 = open("canonicos_alternativos__exon_up", 'w')
	alt_can_7 = open("canonicos_alternativos__exon_down_", 'w') 
	alt_can_8 = open("canonicos_alternativos__exon_down", 'w')		
	
	no_can5_1 = open("no_canonicos5_exon5", "w")
	no_can5_2 = open("no_canonicos5_intron5", "w")
	no_can5_3 = open("no_canonicos5_intron3", "w")	
	no_can5_4 = open("no_canonicos5_exon3", "w")
	no_can5_5 = open("no_canonicos5_intron_up", "w")
	no_can5_6 = open("no_canonicos5_exon_up", "w")
	no_can5_7 = open("no_canonicos5_exon_down", "w")	
	no_can5_8 = open("no_canonicos5_intron_down", "w")	
	
	no_can3_1 = open("no_canonicos3_exon5", "w")
	no_can3_2 = open("no_canonicos3_intron5", "w")
	no_can3_3 = open("no_canonicos3_intron3", "w")	
	no_can3_4 = open("no_canonicos3_exon3", "w")
	no_can3_5 = open("no_canonicos3_intron_up", "w")
	no_can3_6 = open("no_canonicos3_exon_up", "w")
	no_can3_7 = open("no_canonicos3_exon_down", "w")	
	no_can3_8 = open("no_canonicos3_intron_down", "w")	
	
	U12_1 = open("U12_exon5", 'w') 
	U12_2 = open("U12_intron5", 'w')
	U12_3 = open("U12_intron3", 'w') 
	U12_4 = open("U12_exon3", 'w')	

	
	for row in csv.reader(open(final_table), delimiter = '\t'):
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = int(row[5])
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
		non_canonical_filter = row[26].split(",")

		junction_5 = 0
		junction_3 = 0

		if strand == "+":

			if exons_junctions[chr + strand + str(istart)]!= set([]):
				junction_5 = int(max(exons_junctions[chr + strand + str(istart)], key=lambda x:x[1])[0][2])  # selecional el exon con el mayor coverage y luego se indexa para obtener la coordenada genomica

			if exons_junctions[chr + strand + str(iend)]!= set([]):
				junction_3 = int(max(exons_junctions[chr + strand + str(iend)], key=lambda x:x[1])[0][2])

		elif strand == "-":

			if exons_junctions[chr + strand + str(istart)]!= set([]):
				junction_3 = int(max(exons_junctions[chr + strand + str(istart)], key=lambda x:x[1])[0][2])

			if exons_junctions[chr + strand + str(iend)]!= set([]):
				junction_5 = int(max(exons_junctions[chr + strand + str(iend)], key=lambda x:x[1])[0][2])						


		len_up_exon = 1
		len_down_exon =	1			
		len_tag_exon5 = 1
		len_tag_exon3 = 1

		len_tag_intron5 = L
		len_tag_intron3 = L
		len_up_intron = L
		len_down_intron = L
		
		if strand == "+":
			
			while (chr + strand + str(istart - len_tag_exon5) in blocks_starts) == False and len_tag_exon5 <= L:
				len_tag_exon5 += 1				
				
			while (chr + strand + str(iend + len_tag_exon3) in blocks_ends) == False and len_tag_exon3 <= L:
				len_tag_exon3 += 1
				
			while (chr + strand + str(junction_3 - len_down_exon) in blocks_starts) == False and len_down_exon <= L:
				len_down_exon += 1				
				
			while (chr + strand + str(junction_5 + len_up_exon) in blocks_ends) == False and len_up_exon <= L:
				len_up_exon += 1		
		elif strand == "-":
			
			while (chr + strand + str(iend + len_tag_exon5) in blocks_ends) == False and len_tag_exon5 <= L:
				len_tag_exon5 += 1	
			
			while (chr + strand + str(istart - len_tag_exon3) in blocks_starts) == False and len_tag_exon3 <= L:
				len_tag_exon3 += 1

			while (chr + strand + str(junction_5 + len_down_exon) in blocks_ends) == False and len_down_exon <= L:
				len_down_exon += 1	
			
			while (chr + strand + str(junction_3 - len_up_exon) in blocks_starts) == False and len_up_exon <= L:
				len_up_exon += 1
		
		if L > ilength:
			len_tag_intron5 = ilength
			len_tag_intron3 = ilength
		
		intron_up = Genome[chr][junction_5-len_up_intron : junction_5]
		exon_up = Genome[chr][junction_5 : junction_5+len_up_exon]	

		exon5 = Genome[chr][istart-len_tag_exon5:istart]
		intron5 = Genome[chr][istart:istart+len_tag_intron5]
		intron3 = Genome[chr][iend-len_tag_intron3:iend]
		exon3 = Genome[chr][iend:iend+len_tag_exon3]
		
		exon_down = Genome[chr][junction_3-len_down_exon : junction_3]
		intron_down = Genome[chr][junction_3 : junction_3+len_down_intron]	
			
		if strand == "-":
			exon_up = Genome[chr][junction_5-len_down_intron : junction_5].reverse_complement()
			intron_up = Genome[chr][junction_5 : junction_5+len_down_exon].reverse_complement()
			
			exon3 = Genome[chr][istart-len_tag_exon3:istart].reverse_complement()
			intron3 = Genome[chr][istart:istart+len_tag_intron3].reverse_complement()
			intron5 = Genome[chr][iend-len_tag_intron5:iend].reverse_complement()
			exon5 = Genome[chr][iend:iend+len_tag_exon5].reverse_complement()

			intron_down = Genome[chr][junction_3-len_up_exon : junction_3].reverse_complement()
			exon_down = Genome[chr][junction_3 : junction_3+len_up_intron].reverse_complement()
		
		#print intron_up, exon_up, exon_down, intron_down, strand	
		
		intron_up = ">" +  intron + "\n" + str(intron_up).upper()[:L] + "\n" 
		exon_up = ">" +  intron + "\n" + str(exon_up).upper()[:L] + "\n"
		
		exon5 = ">" +  intron + "\n" + str(exon5).upper()[:L] + "\n"
		intron5 = ">" +  intron + "\n" + str(intron5).upper()[:L] + "\n"
		intron3 = ">" +  intron + "\n" + str(intron3).upper()[:L] + "\n"					
		exon3 = ">" +  intron + "\n" + str(exon5).upper()[:L] + "\n"			

		exon_down = ">" +  intron + "\n" + str(exon_down).upper()[:L] + "\n"
		intron_down = ">" +  intron + "\n" + str(intron_down).upper()[:L] + "\n"


		if (dn_type_score >= score_U2_GTAG or dn_type_score >= score_U12_GTAG or dn_type_score >= score_U12_ATAC) and (dn=='GTAG' or dn=='GCAG' or dn=='ATAC') and (ilength >= 40):
				
			if len_tag_exon5 > 8:			
				can_1.write(exon5)
			if len_tag_intron5 > 8:							
				can_2.write(intron5)			
				
			if len_tag_intron3 > 8:
				can_3.write(intron3)
			if len_tag_exon3 > 8:	
				can_4.write(exon3)

			if len_up_intron > 8:
				can_5.write(intron_up)
			if len_up_exon > 8:
				can_6.write(exon_up)
			if len_down_exon > 8:
				can_7.write(exon_down)
			if len_down_intron > 8:
				can_8.write(intron_down)

			
		elif (((dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG and dn[2:]=="AG" and dn!="GCAG") or (dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG and dn[2:]=="AG" and dn!="GCAG") or (dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC and dn[2:]=="AC")) and (dn!='GTAG' or dn!='GCAG' or dn!='ATAC') and (ilength >= 40)):		  

			if len_tag_exon5 > 8:			
				no_can5_1.write(exon5)
			if len_tag_intron5 > 8:							
				no_can5_2.write(intron5)
									
			if len_tag_intron3 > 8:
				no_can5_3.write(intron3)
			if len_tag_exon3 > 8:	
				no_can5_4.write(exon3)

			if len_up_intron > 8:
				no_can5_5.write(intron_up)
			if len_up_exon > 8:
				no_can5_6.write(exon_up)
			if len_down_exon > 8:
				no_can5_7.write(exon_down)
			if len_down_intron > 8:
				no_can5_8.write(intron_down)



		elif (((dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG and dn[:2]=="GT") or (dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG and dn[:2]=="AG" and dn!="GCAG") or (dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC and dn[:2]=="AT"))  and (dn!='GTAG' or dn!='GCAG' or dn!='ATAC') and (ilength >= 40)):

			if len_tag_exon5 > 8:			
				no_can3_1.write(exon5)
			if len_tag_intron5 > 8:							
				no_can3_2.write(intron5)	
				
			if len_tag_intron3 > 8:
				no_can3_3.write(intron3)
			if len_tag_exon3 > 8:	
				no_can3_4.write(exon3)

			if len_up_intron > 8:
				no_can3_5.write(intron_up)
			if len_up_exon > 8:
				no_can3_6.write(exon_up)
			if len_down_exon > 8:
				no_can3_7.write(exon_down)
			if len_down_intron > 8:
				no_can3_8.write(intron_down)

		if (dn_type_score >= score_U2_GTAG or dn_type_score >= score_U12_GTAG and dn_type_score >= score_U12_ATAC) and (dn=='GTAG' or dn=='GCAG' or dn=='ATAC') and (intron_retention_exon != ["NO"] or skipped_exons_names != ["NO"]  or alt_introns != ["NO"]  or alt_no_skipper_introns != ["NO"]  or alt_skipper_introns != ["NO"]  or alt_exon_variant_introns != ["NO"] and (ilength >= 40) ):

			if len_tag_exon5 > 8:			
				alt_can_1.write(exon5)
			if len_tag_intron5 > 8:							
				alt_can_2.write(intron5)
				
			if len_tag_intron3 > 8:
				alt_can_3.write(intron3)
			if len_tag_exon3 > 8:	
				alt_can_4.write(exon3)

			if len_up_intron > 8:
				alt_can_5.write(intron_up)
			if len_up_exon > 8:
				alt_can_6.write(exon_up)
			if len_down_exon > 8:
				alt_can_7.write(exon_down)
			if len_down_intron > 8:
				alt_can_8.write(intron_down)

		if (dn_type=="U12_ATAC" or dn_type =="U12_GTAG") and (dn_type_score >= score_U12_GTAG or dn_type_score >= score_U12_ATAC) and (dn=='GTAG' or dn=='ATAC') and (ilength >= 40):
				
			if len_tag_exon5 > 8:			
				U12_1.write(exon5)
			if len_tag_intron5 > 8:							
				U12_2.write(intron5)			
				
			if len_tag_intron3 > 8:
				U12_3.write(intron3)
			if len_tag_exon3 > 8:	
				U12_4.write(exon3)
							
						
	can_1.close()
	can_2.close()	
	can_3.close()
	can_4.close()
	
	alt_can_1.close()
	alt_can_2.close()	
	alt_can_3.close()
	alt_can_4.close()	

	no_can5_1.close()
	no_can5_2.close()
	no_can5_3.close()
	no_can5_4.close()
	
	no_can3_1.close()
	no_can3_2.close()
	no_can3_3.close()
	no_can3_4.close()



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	Exonextractor(sys.argv[2], sys.argv[3])
	main(sys.argv[3], int(sys.argv[4])) 
