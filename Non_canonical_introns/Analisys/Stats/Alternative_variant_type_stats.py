import sys
import csv


def main(Final_table):

	intron_retention_exon_count = 0
	alt_introns_count = 0
	alt_no_skipper_introns_count = 0
	alt_skipper_introns_count = 0
	alt_exon_variant_introns_count = 0
	shift_count = 0
	non_canonical_shift_count = 0

	for row in csv.reader(open(Final_table), delimiter = '\t'):
		
		gene, intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, tissues_coverage, n_tissues, tissues_names, intron_retention_exon, skipped_exons_names, alt_introns, alt_no_skipper_introns, alt_skipper_introns, alt_exon_variant_introns, shift, non_canonical_shift = row
		
		istart = int(istart)
		iend = int(iend)
		ilength = int(ilength)
		dn_type_score = float(dn_type_score)
		bodymap_coverage = int(bodymap_coverage)
		gm12878_coverage = int(gm12878_coverage)
		hg19_cDNA_coverage = int(hg19_cDNA_coverage)
		hg19_EST_coverage  = int(hg19_EST_coverage)
		mm9_cDNA_coverage = int(mm9_cDNA_coverage)
		mm9_EST_coverage = int(mm9_EST_coverage)
		genecode_coverage = int(genecode_coverage)
		n_tissues = int(n_tissues)
		intron_retention_exon = intron_retention_exon.split(",")
		skipped_exons_names = skipped_exons_names.split(",")
		alt_introns = alt_introns.split(",")
		alt_no_skipper_introns = alt_no_skipper_introns.split(",")
		alt_skipper_introns = alt_skipper_introns.split(",")
		alt_exon_variant_introns = alt_exon_variant_introns.split(",")
		shift = shift.split(",")
		non_canonical_shift = non_canonical_shift.split(",")
		
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
		
			if intron_retention_exon!=['NO']:
				intron_retention = False
				
				for i in intron_retention_exon:
					if i.split("|")[1]!="0.0":
						intron_retention = True
					
					if intron_retention:			
						intron_retention_exon_count += 1

			if alt_introns!=['NO']:
				for i in alt_introns:
					if i.split("|")[2]!="0.0":
						alt_introns_count += 1
		
			if alt_no_skipper_introns!=['NO']:
				for i in alt_no_skipper_introns:
					if i.split("|")[2]!="0.0":
						alt_no_skipper_introns_count += 1		

			if alt_skipper_introns!=['NO']:
				for i in alt_skipper_introns:
					if i.split("|")[2]!="0.0":			
						alt_skipper_introns_count += 1

			if alt_exon_variant_introns!=['NO']:
				for i in alt_exon_variant_introns:
					if i.split("|")[2]!="0.0":				
						alt_exon_variant_introns_count += 1

			if shift!=['NO']:
				for i in shift:
					if i.split("|")[2]!="0.0":		
						shift_count += 1
					
			if non_canonical_shift!=['NO']:
				for i in non_canonical_shift:
					if i.split("|")[2]!="0.0":
						non_canonical_shift_count += 1
		
	print "intron_retention_exon_count, alt_introns_count, alt_no_skipper_introns_count, alt_skipper_introns_count, alt_exon_variant_introns_count, shift_count, non_canonical_shift_count "	
	print intron_retention_exon_count, alt_introns_count, alt_no_skipper_introns_count, alt_skipper_introns_count, alt_exon_variant_introns_count, shift_count, non_canonical_shift_count     
				
			
if __name__ == '__main__':
	main(sys.argv[1])	
