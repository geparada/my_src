import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


Genome = {}
blocks_starts = set([])
blocks_ends = set([])

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
	
	
def Exonextractor (gencode_BED12):
	""" Genera dicionarios de coodenadas de exones que estan presentes en gencode """
		
	reader = csv.reader(open(gencode_BED12), delimiter = '\t')
	
	
	for row in reader:
		chr = row[0]
		start = int(row[1])
		end = row[2]
		name = row[3]
		strand = row[5]
		block_sizes = map (int, row[10].strip(",").split(","))		
		tstarts_0 = map (int, row[11].strip(",").split(",")) 
		
		for t, b in zip(tstarts_0[1:-1], block_sizes[1:-1]):      #No se toman encuenta exones terminales ya que hay muchas isoformas con 5' o 3' alternativos que generan interrupciones de los tags
		
			blocks_starts.add(chr + strand + str(start + t))
			blocks_ends.add(chr + strand + str(start + t + b))



def main (final_table, L):
	""" Extrae secuencias de SJ para analisis de motivos
	La idea no es mezclar secuenencias exonicas con intronicas
	Los tags se interrumpen si es que se encuentran con otro splice juntion, de manera que siempre se escoge el exon mas corto"""
	
	reader = csv.reader(open(final_table), delimiter = '\t')
	
	can_1 = open("canonicos_exon5", 'w') 
	can_2 = open("canonicos_intron5", 'w')
	can_3 = open("canonicos_intron3", 'w') 
	can_4 = open("canonicos_exon3", 'w')
	
	alt_can_1 = open("canonicos_alternativos_exon5", 'w') 
	alt_can_2 = open("canonicos_alternativos_intron5", 'w')
	alt_can_3 = open("canonicos_alternativos_intron3", 'w') 
	alt_can_4 = open("canonicos_alternativos_exon3", 'w')	
	
	no_can5_1 = open("no_canonicos5_exon5", "w")
	no_can5_2 = open("no_canonicos5_intron5", "w")
	no_can5_3 = open("no_canonicos5_intron3", "w")	
	no_can5_4 = open("no_canonicos5_exon3", "w")
	
	no_can3_1 = open("no_canonicos3_exon5", "w")
	no_can3_2 = open("no_canonicos3_intron5", "w")
	no_can3_3 = open("no_canonicos3_intron3", "w")	
	no_can3_4 = open("no_canonicos3_exon3", "w")
	
	U12_1 = open("U12_exon5", 'w') 
	U12_2 = open("U12_intron5", 'w')
	U12_3 = open("U12_intron3", 'w') 
	U12_4 = open("U12_exon3", 'w')	

	
	for row in reader:
		gene, intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, tissues_coverage, n_tissues, tissues_names, intron_retention_exon, skipped_exons_names, alt_introns, alt_no_skipper_introns, alt_skipper_introns, alt_exon_variant_introns, shift, non_canonical_shift = row
		
		istart = int(istart)
		iend = int(iend)
		ilength = int(ilength)
		dn_type_score = float(dn_type_score)
		bodymap_coverage = int(bodymap_coverage)
		gm12878_coverage = int(gm12878_coverage)
		hg19_cDNA_coverage = int(hg19_cDNA_coverage)
		hg19_EST_coverage  = int(hg19_EST_coverage)
		mm9_cDNA_coverage = int(mm9_cDNA_coverage)
		mm9_EST_coverage = int(mm9_EST_coverage)
		genecode_coverage = int(genecode_coverage)
		n_tissues = int(n_tissues)
		intron_retention_exon = intron_retention_exon.split(",")
		skipped_exons_names = skipped_exons_names.split(",")
		alt_introns = alt_introns.split(",")
		alt_no_skipper_introns = alt_no_skipper_introns.split(",")
		alt_skipper_introns = alt_skipper_introns.split(",")
		alt_exon_variant_introns = alt_exon_variant_introns.split(",")
		shift = shift.split(",")
		non_canonical_shift = non_canonical_shift.split(",")	
		
				
		len_tag_exon5 = 1
		len_tag_exon3 = 1
		len_tag_intron5 = L
		len_tag_intron3 = L
		
		if strand == "+":
			
			while (chr + strand + str(istart - len_tag_exon5) in blocks_starts) == False and len_tag_exon5 <= L:
				len_tag_exon5 += 1				
				
			while (chr + strand + str(iend + len_tag_exon3) in blocks_ends) == False and len_tag_exon3 <= L:
				len_tag_exon3 += 1
		
		elif strand == "-":
			
			while (chr + strand + str(iend + len_tag_exon5) in blocks_ends) == False and len_tag_exon5 <= L:
				len_tag_exon5 += 1	
			
			while (chr + strand + str(istart - len_tag_exon3) in blocks_starts) == False and len_tag_exon3 <= L:
				len_tag_exon3 += 1
		
		if L > ilength:
			len_tag_intron5 = ilength
			len_tag_intron3 = ilength

		exon5 = Genome[chr][istart-len_tag_exon5:istart]
		intron5 = Genome[chr][istart:istart+len_tag_intron5]
		intron3 = Genome[chr][iend-len_tag_intron3:iend]
		exon3 = Genome[chr][iend:iend+len_tag_exon3]		
			
		if strand == "-":
			exon3 = Genome[chr][istart-len_tag_exon3:istart].reverse_complement()
			intron3 = Genome[chr][istart:istart+len_tag_intron3].reverse_complement()
			intron5 = Genome[chr][iend-len_tag_intron5:iend].reverse_complement()
			exon5 = Genome[chr][iend:iend+len_tag_exon5].reverse_complement()
		
		exon5 = str(exon5).upper()[:L]
		intron5 = str(intron5).upper()[:L]
		intron3 = str(intron3).upper()[:L]					
		exon3 = str(exon5).upper()[:L]			

		
		exon5 = ">" +  intron + "\n" + str(exon5).upper() + "\n"
		intron5 = ">" +  intron + "\n" + str(intron5).upper() + "\n"
		intron3 = ">" +  intron + "\n" + str(intron3).upper() + "\n"
		exon3 = ">" +  intron + "\n" + str(exon3).upper() + "\n"
		

		donor = ""
		aceptor = ""

		if strand == "+":
			donor = chr + strand + str(istart)  + "SJ5"
			aceptor = chr + strand + str(iend)  + "SJ3" 

		else:
			donor = chr + strand + str(iend) + "SJ5"
			aceptor = chr + strand + str(istart) + "SJ3" 					



		if (dn_type_score >= score_U2_GTAG or dn_type_score >= score_U12_GTAG or dn_type_score >= score_U12_ATAC) and (dn=='GTAG' or dn=='GCAG' or dn=='ATAC') and (ilength >= 40):
				
			if len_tag_exon5 > 8:			
				can_1.write(exon5)
			if len_tag_intron5 > 8:							
				can_2.write(intron5)			
				
			if len_tag_intron3 > 8:
				can_3.write(intron3)
			if len_tag_exon3 > 8:	
				can_4.write(exon3)
			
		elif (((dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG and dn[2:]=="AG" and dn!="GCAG") or (dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG and dn[2:]=="AG" and dn!="GCAG") or (dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC and dn[2:]=="AC")) and (dn!='GTAG' or dn!='GCAG' or dn!='ATAC') and (ilength >= 40)):		  

			if len_tag_exon5 > 8:			
				no_can5_1.write(exon5)
			if len_tag_intron5 > 8:							
				no_can5_2.write(intron5)
									
			if len_tag_intron3 > 8:
				no_can5_3.write(intron3)
			if len_tag_exon3 > 8:	
				no_can5_4.write(exon3)


		elif (((dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG and dn[:2]=="GT") or (dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG and dn[:2]=="AG" and dn!="GCAG") or (dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC and dn[:2]=="AT"))  and (dn!='GTAG' or dn!='GCAG' or dn!='ATAC') and (ilength >= 40)):

			if len_tag_exon5 > 8:			
				no_can3_1.write(exon5)
			if len_tag_intron5 > 8:							
				no_can3_2.write(intron5)	
				
			if len_tag_intron3 > 8:
				no_can3_3.write(intron3)
			if len_tag_exon3 > 8:	
				no_can3_4.write(exon3)

		if (dn_type_score >= score_U2_GTAG or dn_type_score >= score_U12_GTAG and dn_type_score >= score_U12_ATAC) and (dn=='GTAG' or dn=='GCAG' or dn=='ATAC') and (intron_retention_exon != ["NO"] or skipped_exons_names != ["NO"]  or alt_introns != ["NO"]  or alt_no_skipper_introns != ["NO"]  or alt_skipper_introns != ["NO"]  or alt_exon_variant_introns != ["NO"] and (ilength >= 40) ):

			if len_tag_exon5 > 8:			
				alt_can_1.write(exon5)
			if len_tag_intron5 > 8:							
				alt_can_2.write(intron5)
				
			if len_tag_intron3 > 8:
				alt_can_3.write(intron3)
			if len_tag_exon3 > 8:	
				alt_can_4.write(exon3)

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
	Exonextractor(sys.argv[2])
	main(sys.argv[3], int(sys.argv[4])) 
