import sys
import csv

def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')
	
	for row in reader1:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = row[3]
		iend = row[4]
		ilength = row[5]
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
		
#		if row[19]!='NO':
		print "\t".join([chr, istart, iend, intron, str(0), strand])
			

if __name__ == '__main__':
	main(sys.argv[1])
