import sys
import csv

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')

	intron_retention_exon_count = 0
	alt_introns_count = 0
	alt_no_skipper_introns_count = 0
	alt_skipper_introns_count = 0
	alt_exon_variant_introns_count = 0
	shift_count = 0
	non_canonical_shift_count = 0

	constitutives = 0
	alternatives = 0

	for row in reader1:
		
		gene = row[0]
		row = row[1:]
		
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
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
		non_canonical_shift = row[26].split(",")
		
		cons = True
#		print row

#		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
		if dn=="GTAG" or dn=="GCAG" or dn=="ATAC":
			
			
		
			if intron_retention_exon!=['NO']:
				intron_retention = False
				
				for i in intron_retention_exon:
					if i.split("|")[1]!="0.0":
						intron_retention = True
					
					if intron_retention:			
						intron_retention_exon_count += 1
						cons = False


			if alt_introns!=['NO']:
				for i in alt_introns:
					if i.split("|")[2]!="0.0":
						alt_introns_count += 1
						cons = False
		
			if alt_no_skipper_introns!=['NO']:
				for i in alt_no_skipper_introns:
					if i.split("|")[2]!="0.0":
						alt_no_skipper_introns_count += 1
						cons = False

			if alt_skipper_introns!=['NO']:
				for i in alt_skipper_introns:
					if i.split("|")[2]!="0.0":			
						alt_skipper_introns_count += 1
						cons = False

			if alt_exon_variant_introns!=['NO']:
				for i in alt_exon_variant_introns:
					if i.split("|")[2]!="0.0":				
						alt_exon_variant_introns_count += 1
						cons = False

			if shift!=['NO']:
				for i in shift:
					if i.split("|")[2]!="0.0":		
						shift_count += 1
					
			if non_canonical_shift!=['NO']:
				for i in non_canonical_shift:
					if i.split("|")[2]!="0.0":
						non_canonical_shift_count += 1


			if cons:
				 constitutives += 1
			else:
				alternatives += 1


	print "Constitutives", constitutives, percent(constitutives, constitutives + alternatives)
	print "Alternatives",  alternatives, percent(alternatives, constitutives + alternatives)
		
#	print "intron_retention_exon_count, alt_introns_count, alt_no_skipper_introns_count, alt_skipper_introns_count, alt_exon_variant_introns_count, shift_count, non_canonical_shift_count "	
#	print intron_retention_exon_count, alt_introns_count, alt_no_skipper_introns_count, alt_skipper_introns_count, alt_exon_variant_introns_count, shift_count, non_canonical_shift_count     
				
			
if __name__ == '__main__':
	main(sys.argv[1])	
