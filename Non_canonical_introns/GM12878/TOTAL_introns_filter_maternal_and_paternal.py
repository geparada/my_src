import sys
import csv

def main(paternal, maternal):
	reader1 = csv.reader(open(paternal), delimiter = '\t')
	reader2 = csv.reader(open(maternal), delimiter = '\t')
	
#	chr = row[0]		
#	start = row[1]
#	end = row[2]
#	intron = row[3]
#	coverage = row[4]
#	strand = row[5]
#	length = start - end
#	dn = row[8]	
	
	non_can_paternal_list = []
	
	for row in reader1:

		
		chr = row[0]
		start = row[1]
		end = row[2]
		strand = row[5]
		intron = chr + ":" + start + strand + end
		
		coverage_P = row[4]
		dn_P = row[8]
		non_can_paternal_list.append((intron, [dn_P,coverage_P]))
	
	non_can_paternal_dict = dict(non_can_paternal_list)
	
	for row in reader2:
		chr = row[0]
		start = row[1]
		end = row[2]
		strand = row[5]
		intron = chr + ":" + start + strand + end
		length = int(end) - int(start) 
		
		coverage_M = row[4]
		dn_M = row[8]
		
		try:
		
			dn_P = non_can_paternal_dict[intron][0]
			coverage_P = non_can_paternal_dict[intron][1]
		
			
			#if dn_M != 'GTAG' and dn_M != 'GCAG' and dn_M != 'ATAC' and dn_M == dn_P:
			if dn_M == dn_P:			
			
				print intron, min([int(coverage_M),int(coverage_P)]), chr, strand, start, end, length, dn_P
		
		except KeyError:
			pass
		
		
		
	
	
	

	





if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2]) 
