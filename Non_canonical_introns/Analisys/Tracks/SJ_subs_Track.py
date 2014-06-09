import sys
import csv
from collections import defaultdict
import re

tissues_coverage = {}
intron_tissue_coverage = defaultdict(list) # key = intron; contenido = lista de tuplas (nombre del tejido, coverage)

bodymap_out = open('16_tissues_mixture.bed', 'w')
GM12878_out = open('GM12878.bed', 'w')
adipose_out = open('adipose.bed', 'w')
adrenal_out = open('adrenal.bed', 'w')
brain_out = open('brain.bed', 'w')
breast_out = open('breast.bed', 'w')
colon_out = open('colon.bed', 'w')
heart_out = open('heart.bed', 'w')
kidney_out = open('kidney.bed', 'w')
liver_out = open('liver.bed', 'w')
lung_out = open('lung.bed', 'w')
lymph_node_out = open('lymph_node.bed', 'w')
ovary_out = open('ovary.bed', 'w')
prostate_out = open('prostate.bed', 'w')
skeletal_muscle_out = open('skeletal_muscle.bed', 'w')
testes_out = open('testes.bed', 'w')
thyroid_out = open('thyroid.bed', 'w')
white_blood_cells_out = open('white_blood_cells.bed', 'w')


def tissue_reader (tissue, name):
	
	for row in csv.reader(open(tissue), delimiter = ' '):
		intron = row[0]
		coverage = int(row[1])
		intron_tissue_coverage[intron].append((name, coverage))


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
	
	for i in intron_tissue_coverage.items():
		intron = i[0]
		tissues_coverage_tuples = i[1]
		tissues_coverage[intron] = dict(tissues_coverage_tuples)
	
	
	

def track_gen (intron, dn, dn_type, dn_type_score, RGB, coverage, out_file):
	""" Genera los BED12 """

	if coverage > 0:
		
		chr, istart, iend = re.findall(r"[\w']+", intron)
		istart = int(istart)
		iend = int(iend)
		
		strand = ""
		if "+" in intron:
			strand = "+"
		if "-" in intron:
			strand = "-"

		score_U2_GTAG = float(63.2891027914)
		score_U12_ATAC = float(60.9280810964)
		score_U12_GTAG = float(61.4553595446)
					
		ID = dn[:2] + "-" + dn[2:]  + "[" + str(coverage) + "]"
		
		start = istart-8
		end = iend+8

		blockCount = "2"
		blockSizes = "8,8"
		blockStarts = "0" + "," + str(iend-start)
			
		start = istart-8
		end = iend+8
			
		if dn!="GTAG" and dn!="ATAC" and dn!="GCAG":
			
			red = 0
			green = 0 
			blue = 0
				
			if dn_type == "U2_GTAG":
				if dn_type_score >= score_U2_GTAG :
					green = int(100 + ((dn_type_score-score_U2_GTAG)*155)/(100-score_U2_GTAG))
				else:
					red = int(100 + ((score_U2_GTAG-dn_type_score)*155)/(score_U2_GTAG-20))

					
			elif dn_type == "U12_ATAC":
				if dn_type_score >= score_U12_ATAC:
					green = int(100 + ((dn_type_score-score_U2_GTAG)*155)/(100-score_U2_GTAG))
				else:
					red = int(100 + ((score_U12_ATAC-dn_type_score)*155)/(score_U12_ATAC-20))
			
					
			elif dn_type == "U12_GTAG":
				if dn_type_score >= score_U12_GTAG:
					green = int(100 + ((dn_type_score-score_U2_GTAG)*155)/(100-score_U2_GTAG))
				else:
					red = int(100 +  ((score_U12_GTAG-dn_type_score)*155)/(score_U12_GTAG-20))
			
			RGB = ",".join([str(red),str(green),str(blue)]) 
			
		BED = [chr, str(start), str(end ), ID, "0", strand, str(start), str(end), RGB, blockCount, blockSizes, blockStarts]
		out_file.write("\t".join(BED)+"\n") 




def main(Final_table):	

	for row in csv.reader(open(Final_table), delimiter = '\t'):  
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = float(row[8])
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		TOTAL_tissues_coverage = int(row[16])
		n_tissues = row[17]
		tissues = row[18]

		adipose_coverage = 0
		adrenal_coverage = 0
		brain_coverage = 0
		breast_coverage = 0
		colon_coverage = 0
		heart_coverage = 0
		kidney_coverage = 0
		liver_coverage = 0
		lung_coverage = 0
		lymph_node_coverage = 0
		ovary_coverage = 0
		prostate_coverage = 0
		skeletal_muscle_coverage = 0
		testes_coverage = 0
		thyroid_coverage = 0
		white_blood_cells_coverage = 0
				
		try:		
			adipose_coverage = tissues_coverage[intron]["adipose"]
		except KeyError:
			pass	

		try:		
			adrenal_coverage = tissues_coverage[intron]["adrenal"]
		except KeyError:
			pass		

		try:		
			brain_coverage = tissues_coverage[intron]["brain"]
		except KeyError:
			pass		

		try:		
			breast_coverage = tissues_coverage[intron]["breast"]
		except KeyError:
			pass

		try:		
			colon_coverage = tissues_coverage[intron]["colon"]
		except KeyError:
			pass
			
		try:		
			heart_coverage = tissues_coverage[intron]["heart"]
		except KeyError:
			pass
			
		try:		
			kidney_coverage = tissues_coverage[intron]["kidney"]
		except KeyError:
			pass

		try:		
			liver_coverage = tissues_coverage[intron]["liver"]
		except KeyError:
			pass

		try:		
			lung_coverage = tissues_coverage[intron]["lung"]
		except KeyError:
			pass	

		try:		
			lymph_node_coverage = tissues_coverage[intron]["lymph_node"]
		except KeyError:
			pass		

		try:		
			ovary_coverage = tissues_coverage[intron]["ovary"]
		except KeyError:
			pass		

		try:		
			prostate_coverage = tissues_coverage[intron]["prostate"]
		except KeyError:
			pass

		try:		
			skeletal_muscle_coverage = tissues_coverage[intron]["skeletal_muscle"]
		except KeyError:
			pass
			
		try:		
			testes_coverage = tissues_coverage[intron]["testes"]
		except KeyError:
			pass
			
		try:		
			thyroid_coverage = tissues_coverage[intron]["thyroid"]
		except KeyError:
			pass

		try:		
			white_blood_cells_coverage = tissues_coverage[intron]["white_blood_cells"]
		except KeyError:
			pass
		
		
		track_gen (intron, dn, dn_type, dn_type_score,"142,142,56", bodymap_coverage, bodymap_out)  #olivedrab
		track_gen (intron, dn, dn_type, dn_type_score, "238,130,98", gm12878_coverage, GM12878_out)	#salmon2	
		track_gen (intron, dn, dn_type, dn_type_score, "205,198,115" , adipose_coverage, adipose_out) #khaki3		
		track_gen (intron, dn, dn_type, dn_type_score, "556,142,142", adrenal_coverage, adrenal_out) #sgi teal 	
		track_gen (intron, dn, dn_type, dn_type_score, "255,153,18", brain_coverage, brain_out) #cadmiumyellow
		track_gen (intron, dn, dn_type, dn_type_score, "255,62,150", breast_coverage, breast_out) #violetred1		
		track_gen (intron, dn, dn_type, dn_type_score, "0,255,127", colon_coverage, colon_out)	#springgreen	
		track_gen (intron, dn, dn_type, dn_type_score,"139,10,80", heart_coverage, heart_out) #deeppink4
		track_gen (intron, dn, dn_type, dn_type_score, "192,255,62", kidney_coverage, kidney_out) #olivedrab
		track_gen (intron, dn, dn_type, dn_type_score, "113,113,198", liver_coverage, liver_out) #sgi slateblue		
		track_gen (intron, dn, dn_type, dn_type_score, "0,255,255", lung_coverage, lung_out) #cyan/aqua		
		track_gen (intron, dn, dn_type, dn_type_score, "139,76,57", lymph_node_coverage, lymph_node_out) #salmon4	
		track_gen (intron, dn, dn_type, dn_type_score, "255,131,250", ovary_coverage, ovary_out) #ochid1
		track_gen (intron, dn, dn_type, dn_type_score, "255,193,37", prostate_coverage, prostate_out)	#goldenrod1	
		track_gen (intron, dn, dn_type, dn_type_score, "198,133,133", skeletal_muscle_coverage, skeletal_muscle_out) #sgi salmon		
		track_gen (intron, dn, dn_type, dn_type_score, "16,78,139", testes_coverage, testes_out)	#dodgerblue4													
		track_gen (intron, dn, dn_type, dn_type_score, "255,182,193", thyroid_coverage, thyroid_out) #ligtpink		
		track_gen (intron, dn, dn_type, dn_type_score, "139,131,134", white_blood_cells_coverage, white_blood_cells_out) #lavenderblush4	

#ver cloford.com/resources/colours/500col.html


# python ~/my_src/Analisys/Tracks/SJ_subs_Track.py FINAL_TABLE.alternative_splicing.fold.highfold adipose.TOTAL_SJ_coverage adrenal.TOTAL_SJ_coverage brain.TOTAL_SJ_coverage breast.TOTAL_SJ_coverage colon.TOTAL_SJ_coverage heart.TOTAL_SJ_coverage kidney.TOTAL_SJ_coverage liver.TOTAL_SJ_coverage lung.TOTAL_SJ_coverage lymph_node.TOTAL_SJ_coverage ovary.TOTAL_SJ_coverage prostate.TOTAL_SJ_coverage skeletal_muscle.TOTAL_SJ_coverage testes.TOTAL_SJ_coverage thyroid.TOTAL_SJ_coverage white_blood_cells.TOTAL_SJ_coverage

if __name__ == '__main__':
	tissue_procesing(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17])	
	main(sys.argv[1])



bodymap_out.close()
GM12878_out.close()
adipose_out.close()
adrenal_out.close()
brain_out.close()
breast_out.close()
colon_out.close()
heart_out.close()
kidney_out.close()
liver_out.close()
lung_out.close()
lymph_node_out.close()
ovary_out.close()
prostate_out.close()
skeletal_muscle_out.close()
testes_out.close()
thyroid_out.close()
white_blood_cells_out.close()
