import sys
import csv
from collections import defaultdict


#def main(TOTAL_alineados_al_genoma, TOTAL_pseudogenes, GM12878_SJ, HeLaS3_SJ, HepG2_SJ, HUVEC_SJ, H1hesc_SJ):
def main(TOTAL_alineados_al_genoma, TOTAL_pseudogenes, polyA, RNA_total):	
	
	reads_SJ_genoma = []
	
	reader1 = csv.reader(open(TOTAL_alineados_al_genoma), delimiter = '\t')
	reader2 = csv.reader(open(TOTAL_pseudogenes), delimiter = '\t')
	
	reader3 = csv.reader(open(polyA), delimiter = '\t')
	reader4 = csv.reader(open(RNA_total), delimiter = '\t')

		
	reads_mismatches_genoma = {}
	reads_pseudogenes = defaultdict(list)
	reads_mismatches_pseudogenes = {}
	
	for row in reader1:
		read = row[0]
		flag = int(row[1])
		mismatchs = 0
	
		if flag !=4:
			mismatches = int(row[13].strip("NM:i:"))
			reads_mismatches_genoma[read] = mismatches

	for row in reader2:
		read = row[0]
		flag = int(row[1])
		pseudogen = row[2]
		if flag !=4:
			mismatches = int(row[13].strip("NM:i:"))
			reads_mismatches_pseudogenes[read] = mismatches
			reads_pseudogenes[read].append(pseudogen)
	
		
	print "read", "RNA-seq_library", "SJ", "mismatches_SJ", "mismatches_genoma", "mismatches_pseudogen", "peseudogenes", "seq", "qual"
	
	for row in reader3:
		read = row[0]
		SJ = row[2]		
		mismatches_SJ = int(row[13].strip("NM:i:"))
		mismatches_genoma = "NO"
		mismatches_pseudogen = "NO"
		pseudogen = "NO"
		cell_line = "polyA"
		seq = row[9]
		qual = row[10]
		
		if reads_mismatches_genoma.has_key(read):
			mismatches_genoma = reads_mismatches_genoma[read]
		
		if reads_pseudogenes.has_key(read):
			pseudogen = reads_pseudogenes[read]
			mismatches_pseudogen = reads_mismatches_pseudogenes[read]
		
		if pseudogen!="NO":
			peseudogenes = ",".join(pseudogen)
		else:
			peseudogenes = pseudogen
		
		print read, cell_line, SJ, mismatches_SJ, mismatches_genoma, mismatches_pseudogen, peseudogenes, seq, qual
		
		
	for row in reader4:
		read = row[0]
		SJ = row[2]			
		mismatches_SJ = int(row[13].strip("NM:i:"))
		mismatches_genoma = "NO"
		mismatches_pseudogen = "NO"	
		pseudogen = "NO"
		cell_line = "RNA_total"
		seq = row[9]
		qual = row[10]
		
		if reads_mismatches_genoma.has_key(read):
			mismatches_genoma = reads_mismatches_genoma[read]
		
		if reads_pseudogenes.has_key(read):
			pseudogen = reads_pseudogenes[read]
			mismatches_pseudogen = reads_mismatches_pseudogenes[read]
		
		if pseudogen!="NO":
			peseudogenes = ",".join(pseudogen)
		else:
			peseudogenes = pseudogen
		
		print read, cell_line, SJ, mismatches_SJ, mismatches_genoma, mismatches_pseudogen, peseudogenes, seq, qual






#python ~/my_src/Roberto/Mirror_ENCODE.py TOTAL_alineados_al_genoma.sam TOTAL_pseudogenes.sam Reads_mirror_GM12878_SJs.tabular Reads_mirror_HeLaS3_SJs.tabular Reads_mirror_HepG2_SJs.tabular Reads_mirror_HUVEC_SJs.tabular Reads_mirrors_H1hesc_SJs.tabular



if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2], sys.argv[3],sys.argv[4])
