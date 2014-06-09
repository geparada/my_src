import sys
import csv

def main (non_canonical_final_table):
	""" Genera un BED de las bases de los dinucleotidos """
	
	for row in csv.reader(open(non_canonical_final_table), delimiter = '\t'):
		
		name = row[0]  	
		intron = row[1]
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		ilength = row[6]
		dn = row[7]
		dn_type = row[8]
		dn_type_score = float(row[9])
		bodymap_coverage = int(row[10])
		gm12878_coverage = int(row[11])
		hg19_cDNA_coverage = int(row[12])
		hg19_EST_coverage  = int(row[13])
		mm9_cDNA_coverage = int(row[14])
		mm9_EST_coverage = int(row[15])
		genecode_coverage = int(row[16])
		TOTAL_tissues_coverage = int(row[17])
		n_tissues = row[18]	
		intron_retention_exon = row[19].split(",")
		skipped_exons_names = row[20].split(",")
		alt_introns = row[21].split(",")
		alt_no_skipper_introns = row[22].split(",")
		alt_skipper_introns = row[23].split(",")
		alt_exon_variant_introns = row[24].split(",")
		shift = row[25].split(",")
		non_canonical_shift = row[26].split(",")
		
		if strand == "+":
		
			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, istart, istart+1, intron + "_1", 0, strand) 
			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, istart+1, istart+2, intron + "_2", 0, strand)
			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, iend-2, iend-1, intron + "_4", 0, strand)					
			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, iend-1, iend, intron + "_3", 0, strand)
		
		if strand == "-":

			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, iend-1, iend, intron + "_3", 0, strand)
			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, iend-2, iend-1, intron + "_4", 0, strand)		
			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, istart+1, istart+2, intron + "_2", 0, strand)		
			print '%s\t%s\t%s\t%s\t%s\t%s' %(chr, istart, istart+1, intron + "_1", 0, strand) 

			
					


if __name__ == '__main__':
	main(sys.argv[1])
