import sys
import csv

def main(paternal, maternal):
	csv.field_size_limit(1000000000)
	reader1 = csv.reader(open(paternal), delimiter = ' ')
	reader2 = csv.reader(open(maternal), delimiter = ' ')
	
	intron_dict = {}
	
	for row in reader1:
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
		reads = set(row[10].split(","))
		
		intron_dict[intron] = [coverage, reads]
	
	for row in reader2:
		intron = row[0]
		maternal_coverage = int(row[1])
		chr = row[2]
		strand = row[3]
		istart = row[4]
		iend = row[5]
		ilength = int(row[6])
		dn = row[7]
		dn_type = row[8]
		dn_type_score = row[9]
		maternal_reads = set(row[10].split(","))
		
		try:
			paternal_coverage = intron_dict[intron][0]
			
			if dn != "GTAG" and dn != "GCAG" and dn != "ATAC":

				paternal_reads = intron_dict[intron][1]
				common_reads = maternal_reads & paternal_reads
				
				#if len(common_reads) >= 3:
				print intron, len(common_reads) , chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, ",".join(common_reads)
		
			else:
				print  intron, min(maternal_coverage, paternal_coverage)  , chr, strand, istart, iend, ilength,  dn, dn_type, dn_type_score, 0
		
		except KeyError:
			pass


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
