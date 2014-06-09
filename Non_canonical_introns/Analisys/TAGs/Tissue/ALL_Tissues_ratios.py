import sys
import csv
from collections import defaultdict


ALL = defaultdict(list)

def tissiue_reader (tissiue):
	""" Lee y procesa los ratios de los distintos tejidos """
	
	reader1 = csv.reader(open(tissiue), delimiter = ' ')
	
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
		tissue_coverage = int(row[17])
		bodymap_seq = row[18]
		gm12878_seq = row[19]
		DR = row[20]
		intron_retention_exon = "none" #row[21]
		skipped_exons_names = "none" #row[22]
		alt_introns = row[23].split(",")
		alt_no_skipper_introns = row[24].split(",")
		alt_skipper_introns = row[25].split(",")
		alt_exon_variant_introns = row[26].split(",")
		shift = row[27].split(",")
		non_canonical_filter = row[28].split(",")
		
		add_INC_IC (gene, intron, dn, alt_no_skipper_introns, "alt_no_skipper_intron")
		add_INC_IC (gene, intron, dn, alt_skipper_introns, "alt_skipper_intron")
		add_INC_IC (gene, intron, dn, alt_exon_variant_introns, "alt_exon_variant_intron")
		add_INC_IC (gene, intron, dn, shift, "shift")
		add_INC_IC (gene, intron, dn, non_canonical_filter, "non_canonical_shift")		
	
	

def add_INC_IC (INC_gene, INC_intron, INC_dn, lista_alt, alt_type):
	""" Agrega las parejas de intrones no-canonicos/canonicos al dicionario"""

	for i in lista_alt:
		if i !="NO":
			alt_intron = i.split("|")[0]
			alt_dn = i.split("|")[1]
			ratio_can = i.split("|")[2].split("/")[0]
			ratio_non_can = i.split("|")[2].split("/")[1]
			ratio = ratio_non_can + "/" + ratio_can
			alt_donor = i.split("|")[3]
			alt_aceptor = i.split("|")[4]

			key = INC_gene + '|' + INC_intron + '|' + INC_dn + '|' + alt_type + '|' + alt_donor + '|' + alt_aceptor + '|'  + alt_intron + '|' + alt_dn 
			ALL[key].append(ratio)
	


def main (adipose, adrenal, brain, breast, colon, heart, kidney, liver, lung, lymph_node, ovary, prostate, skeletal_muscle, testes, thyroid,  white_blood_cells):
	""" Genera tabla resumen de todos los ratios entre intrones no canonicos y canonicos en cada uno de los tejidos """

	tissiue_reader(adipose)
	tissiue_reader(adrenal)
	tissiue_reader(brain)
	tissiue_reader(breast)	
	tissiue_reader(colon)
	tissiue_reader(heart)
	tissiue_reader(kidney)
	tissiue_reader(liver)	
	tissiue_reader(lung)	
	tissiue_reader(lymph_node)	
	tissiue_reader(ovary)	
	tissiue_reader(prostate)	
	tissiue_reader(skeletal_muscle)	
	tissiue_reader(testes)
	tissiue_reader(thyroid)	
	tissiue_reader(white_blood_cells)	

	print "gene", "INC", "INC_dn", "alt_type", "alt_donor",  "alt_aceptor", "alt_intron", "alt_dn", "adipose", "adrenal", "brain", "breast", "colon", "heart", "kidney", "liver", "lung", "lymph_node", "ovary", "prostate", "skeletal_muscle", "testes", "thyroid",  "white_blood_cells"  			

	for row in ALL.items():
			
		print row[0].replace("|", " "), " ".join(row[1])
	

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16])
