import sys
import csv
from decimal import *
import math
getcontext().prec = 20
from collections import defaultdict


def main(Final_table):

	black_list = set([])

	for row in csv.reader(open(Final_table), delimiter = '\t'):

		intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, tissues_coverage, n_tissues, tissues_names, intron_retention_exon, skipped_exons_names, alt_introns, alt_no_skipper_introns, alt_skipper_introns, alt_exon_variant_introns, shift, non_canonical_shift = row
		
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
		
		folds = []
		
		if intron_retention_exon!=['NO']:
			for i in intron_retention_exon:
				folds.append( Decimal(i.split("|")[1]))

		if skipped_exons_names!=['NO']:
			for i in skipped_exons_names:
				folds.append( Decimal(i.split("|")[1]))

		if alt_introns!=['NO']:
			for i in alt_introns:
				folds.append( Decimal(i.split("|")[2]))

		if shift!=['NO']:
			for i in shift:
				folds.append( Decimal(i.split("|")[2]))
				
		max_fold = 0
		
		if len(folds)!=0:
			max_fold = max(folds)

		if max_fold >= 20:
			black_list.add(intron)
		
		
	for row in csv.reader(open(Final_table), delimiter = '\t'):
		intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, tissues_coverage, n_tissues, tissues_names, intron_retention_exon, skipped_exons_names, alt_introns, alt_no_skipper_introns, alt_skipper_introns, alt_exon_variant_introns, shift, non_canonical_shift = row
		
		alt_introns = set(alt_introns.split(",")) - black_list
		alt_no_skipper_introns = set(alt_no_skipper_introns.split(",")) - black_list
		alt_skipper_introns = set(alt_skipper_introns.split(",")) - black_list
		alt_exon_variant_introns = set(alt_exon_variant_introns.split(",")) - black_list
		shift = set(shift.split(",")) - black_list
		non_canonical_shift = set(non_canonical_shift.split(",")) - black_list

		if (intron in black_list) == False:
			
			alt_table = [ ','.join(alt_introns), ','.join(alt_no_skipper_introns), ','.join(alt_skipper_introns),  ','.join(alt_exon_variant_introns), ','.join(shift), ','.join(non_canonical_shift)]

			print "\t".join(row[:21]) + "\t" + "\t".join(["NO" if len(x)==0 else x for x in alt_table])


			
			
if __name__ == '__main__':
	main(sys.argv[1])
