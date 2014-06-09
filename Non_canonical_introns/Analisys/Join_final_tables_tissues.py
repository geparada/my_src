import sys
import csv
from collections import defaultdict

total_introns = set([])
tissues = defaultdict(set)
tissues_reads = defaultdict(int)

def dic_gen(file):
	
	dict_info = {}
	
	csv.field_size_limit(1000000000)
	reader = csv.reader(open(file), delimiter = ' ')

	for row in reader:
		
		intron = row[0]
		coverage = int(row[1])
		chr = row[2]
		strand = row[3]
		istart = row[4]
		iend = row[5]
		ilength = int(row[6])
		dn = row[7]
		dn_type = row[8]
		dn_type_score = row[9]
		reads = row[10]

		total_introns.add((intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score))
		dict_info[intron] =  [coverage, ilength]       #Para los cDNA_EST el ilenght es el coverage de EST 		
	
	return dict_info

def search(intron, dict_info):
	try:
		info = dict_info[intron]
	except KeyError:
		info = [0, 0, 0, 0, 0, 0, 0, 0]
	return info
	

def tissue_reader (tissue, name):
	for row in csv.reader(open(tissue), delimiter = ' '):
		intron = row[0]
		coverage = int(row[1])
		tissues_reads[intron] += coverage
		if coverage >= 3:
			tissues[intron].add(name)

def tissue_procesing (adipose, adrenal, brain, breast, colon, heart, kidney, liver, lung, lymph_node, ovary, prostate, skeletal_muscle, testes, thyroid,  white_blood_cells):
	
	tissue_reader(adipose,'adipose')
	tissue_reader(adrenal, 'adrenal')
	tissue_reader(brain, 'brain')
	tissue_reader(breast, 'breast')	
	tissue_reader(colon, 'colon')
	tissue_reader(heart, 'heart')
	tissue_reader(kidney, 'kidney')
	tissue_reader(liver, 'liver')	
	tissue_reader(lung, 'lung')	
	tissue_reader(lymph_node, 'lymph_node')	
	tissue_reader(ovary, 'ovary')	
	tissue_reader(prostate, 'prostate')	
	tissue_reader(skeletal_muscle, 'skeletal_muscle')	
	tissue_reader(testes, 'testes')
	tissue_reader(thyroid, 'thyroid')	
	tissue_reader(white_blood_cells, 'white_blood_cells')			


def main(bodymap, gm12878, hg19_cDNA, hg19_EST, mm9_cDNA, mm9_EST, genecode):
	
	dict_bodymap = 	dic_gen(bodymap)
	dict_gm12878 = 	dic_gen(gm12878)
	dict_hg19_cDNA = dic_gen(hg19_cDNA)
	dict_hg19_EST = dic_gen(hg19_EST)
	dict_mm9_cDNA = dic_gen(mm9_cDNA)
	dict_mm9_EST = dic_gen(mm9_EST)
	dict_genecode = dic_gen(genecode)

	
	for row in total_introns:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = row[3]
		iend = row[4]
		ilength =row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = row[8]
		tissues_coverage = tissues_reads[intron]
		n_tissues = len(tissues[intron])
					
		bodymap_info = search(intron, dict_bodymap)
		gm12878_info = search(intron, dict_gm12878)
		hg19_cDNA_info = search(intron, dict_hg19_cDNA)
		hg19_EST_info = search(intron, dict_hg19_EST)
		mm9_cDNA_info = search(intron, dict_mm9_cDNA)
		mm9_EST_info = search(intron, dict_mm9_EST)
		genecode_info = search(intron, dict_genecode)
		
		bodymap_coverage =  bodymap_info[0]
		gm12878_coverage =	gm12878_info[0]
		hg19_cDNA_coverage = hg19_cDNA_info[0]
		hg19_EST_coverage = hg19_EST_info[0]
		mm9_cDNA_coverage = mm9_cDNA_info[0]
		mm9_EST_coverage = mm9_EST_info[0]
		genecode_coverage =	genecode_info[0]	

		hg_coverage = hg19_cDNA_coverage + hg19_EST_coverage
		mm9_coverage = mm9_cDNA_coverage + mm9_EST_coverage
		tissues_names = ','.join(tissues[intron])
		if tissues_names == "":
			tissues_names = 0
				
#		if (dn!='GTAG' and dn!='GCAG' and dn!='ATAC'): # and (bodymap_coverage>=3 or gm12878_coverage>=3 or hg_coverage>=3):
		
#		if intron == "chr12:15103657-15114470":
		
		if ((bodymap_coverage>=3 and gm12878_coverage>=3) or (bodymap_coverage>=3 and hg_coverage>=3) or (gm12878_coverage>=3 and hg_coverage>=3) or (gm12878_coverage>=3 and n_tissues>=1) or (hg_coverage>=3 and n_tissues>=1) or (len(tissues[intron] - set(["kidney", "thyroid"])) >= 1 and n_tissues>=2) or (len(tissues[intron] - set(["heart", "skeletal_muscle"])) >= 1 and n_tissues>=2)):
		
#		if  (dn!='GTAG' and dn!='GCAG' and dn!='ATAC') and ( ( (bodymap_coverage>=3 and gm12878_coverage>=3) or (bodymap_coverage>=3 and hg_coverage>=3) or (gm12878_coverage>=3 and hg_coverage>=3) ) or ( (bodymap_coverage>=3 and mm9_coverage>=3) or (gm12878_coverage>=3 and mm9_coverage>=3)  ) or (n_tissues >=3) ):   #(intreseccion) or (conservados y RNA-seq) or (en mas de dos tejidos)
			   		
			print intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, tissues_coverage, n_tissues, tissues_names 
	
#kidney, thyroid
#heart, skeletal_muscle 

# python ~/my_src/Analisys/Join_final_tables_tissues.py ../hg19/ALL/introns.final_table.hg19.fixed.deletions.SNPs.tags ../GM12878/NA12878_Joel_Rozowsky/STRANDED/TOTAL/introns.final_table.hg19.fixed.tags.patch_U12_GTAG ../hg19/cDNAS_ESTs/introns.mrna.fixed.deletions.SNPs ../hg19/cDNAS_ESTs/introns.est.fixed.deletions.SNPs ../mm9/New/introns.mrna.fixed.hg19.fixed ../mm9/New/introns.est.fixed.hg19.fixed  ../genecode11/introns.genecode.fixed ../Tissue/tags/adipose.TOTAL_SJ_coverage ../Tissue/tags/adrenal.TOTAL_SJ_coverage ../Tissue/tags/brain.TOTAL_SJ_coverage ../Tissue/tags/breast.TOTAL_SJ_coverage ../Tissue/tags/colon.TOTAL_SJ_coverage ../Tissue/tags/heart.TOTAL_SJ_coverage ../Tissue/tags/kidney.TOTAL_SJ_coverage ../Tissue/tags/liver.TOTAL_SJ_coverage ../Tissue/tags/lung.TOTAL_SJ_coverage ../Tissue/tags/lymph_node.TOTAL_SJ_coverage ../Tissue/tags/ovary.TOTAL_SJ_coverage ../Tissue/tags/prostate.TOTAL_SJ_coverage ../Tissue/tags/skeletal_muscle.TOTAL_SJ_coverage ../Tissue/tags/testes.TOTAL_SJ_coverage ../Tissue/tags/thyroid.TOTAL_SJ_coverage ../Tissue/tags/white_blood_cells.TOTAL_SJ_coverage

if __name__ == '__main__':
	tissue_procesing(sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23])	
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
