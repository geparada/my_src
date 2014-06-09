import sys
import csv

total_introns = set([])

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
				
		if (bodymap_coverage>=3 or gm12878_coverage>=3 or hg_coverage>=3):			
			   		
			print intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage
	




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
