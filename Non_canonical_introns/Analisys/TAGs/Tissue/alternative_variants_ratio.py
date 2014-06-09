import sys
import csv
csv.field_size_limit(1000000000)
from decimal import *
import math
getcontext().prec = 20

#Usar FINAL_TABLE.alternative_splicing.fold.highfold.gene_names




def new_ratios (non_can_coverage, alt_list, coverage_dict, new_alt_list):
	""" Recalcula los nuevos ratios """

	for i in alt_list:
		
		if i !="NO":
			alt_intron = i.split("|")[0]
			alt_dn = i.split("|")[1]
			alt_donor = i.split("|")[3]
			alt_aceptor = i.split("|")[4]
									
			alt_coverage = 0     #esto es lo que se recalcula para cada tejido
				
			try:
				alt_coverage = coverage_dict[alt_intron]
			except KeyError:
				pass
				
			ratio = 0
			
			ratio = str(alt_coverage) + "/" + str(non_can_coverage)
				
#			if non_can_coverage != 0:
#				ratio = round (Decimal(alt_coverage) / Decimal(non_can_coverage), 3)
#			elif alt_coverage != 0:
#				ratio = str(alt_coverage) + "/0"
									
			new_alt_list.append(alt_intron + "|" + alt_dn + "|" + str(ratio) + "|" + alt_donor + "|" + alt_aceptor)
			
			
def list_join (X_list):
	""" Junta los elementos de una lista """
	
	if X_list != []:
		return ",".join(X_list)
	else:
		return "NO"
	
	

def main (alternative_splicing, SJ_coverage):
	"""Calcular el ratio entre el coverage de los intrones no canonicos y sus distitas variantes de splicing"""
	
	coverage_introns = {}
	

	reader1 = csv.reader(open(alternative_splicing), delimiter = '\t')
	reader2 = csv.reader(open(SJ_coverage), delimiter = ' ')
	
	for row in reader2:
		intron = row[0]
		coverage = int(row[1])
		coverage_introns[intron] = coverage		

	for row in reader1:
		
		gene = row[0]
		intron = row[1]
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		ilength = row[6]
		dn = row[7]
		dn_type = row[8]
		dn_type_score = row[9]
		bodymap_coverage = int(row[10])
		gm12878_coverage = int(row[11])
		hg19_cDNA_coverage = int(row[12])
		hg19_EST_coverage  = int(row[13])
		mm9_cDNA_coverage = int(row[14])
		mm9_EST_coverage = int(row[15])
		genecode_coverage = int(row[16])
		bodymap_seq = row[17]
		gm12878_seq = row[18]
		DR = row[19]
		intron_retention_exon = "none" #row[19]
		skipped_exons_names = "none" #row[20]
		alt_introns = row[22].split(",")
		alt_no_skipper_introns = row[23].split(",")
		alt_skipper_introns = row[24].split(",")
		alt_exon_variant_introns = row[25].split(",")
		shift = row[26].split(",")
		non_canonical_filter = row[27].split(",")
		
		new_tissue_coverage = 0
		
		try:
			new_tissue_coverage =  coverage_introns[intron]
		
		except KeyError:
			pass
		
		
		new_alt_introns = []
		new_alt_no_skipper_introns = []
		new_alt_skipper_introns = []
		new_alt_exon_variant_introns = []
		new_shift = []
		new_non_canonical_filter = []
		
		
		new_ratios(new_tissue_coverage, alt_introns, coverage_introns, new_alt_introns)
		new_ratios(new_tissue_coverage, alt_no_skipper_introns, coverage_introns, new_alt_no_skipper_introns)
		new_ratios(new_tissue_coverage, alt_skipper_introns, coverage_introns, new_alt_skipper_introns)
		new_ratios(new_tissue_coverage, alt_exon_variant_introns, coverage_introns, new_alt_exon_variant_introns)
		new_ratios(new_tissue_coverage, shift, coverage_introns, new_shift)
		new_ratios(new_tissue_coverage, non_canonical_filter, coverage_introns, new_non_canonical_filter)
		
			
		print gene, intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage, mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, new_tissue_coverage, bodymap_seq, gm12878_seq, DR, intron_retention_exon, skipped_exons_names, list_join(new_alt_introns), list_join(new_alt_no_skipper_introns), list_join(new_alt_skipper_introns), list_join(new_alt_exon_variant_introns), list_join(new_shift), list_join(new_non_canonical_filter)


										

#ZeroDivisionError:			   
		
		#print intron, dn, dn_type, new_tissue_coverage, intron_retention_exon

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
