import sys
import csv
from collections import defaultdict




def main(Final_table):   #La completa, con canonicos y no canonicos

	U2_U12_like = open(str(sys.argv[1]) + ".U2_U12_like", 'w')
	indef = open(str(sys.argv[1]) + ".indef", 'w')
	indef_shift = open(str(sys.argv[1]) + ".indef.shift", 'w')
	indef_non_shift = open(str(sys.argv[1]) + ".indef.non_shift", 'w')

	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)
	

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


		if dn!="GTAG" and dn!="ATAC" and dn!="GCAG":

			min_ilength = 85
			high_score = 70

			rescue = alt_introns != ['NO'] or (mm9_EST_coverage + mm9_cDNA_coverage) >= 3

			if ilength >= min_ilength:		
		
				if dn_type == "U2_GTAG" and dn_type_score >= high_score:
					U2_U12_like.write("\t".join(row) + "\n")


				elif dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG and rescue:
					U2_U12_like.write("\t".join(row)+ "\n")
					
				elif dn_type == "U12_ATAC" and dn_type_score >= high_score:
					U2_U12_like.write("\t".join(row)+ "\n")

				elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC and rescue:
					U2_U12_like.write("\t".join(row)+ "\n")

					
				elif dn_type == "U12_GTAG" and dn_type_score >= high_score:
					U2_U12_like.write("\t".join(row)+ "\n")

				elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG and rescue:
					U2_U12_like.write("\t".join(row)+ "\n")
					
				else:
					indef.write("\t".join(row)+ "\n")

					if shift == ['NO']:
						indef_non_shift.write("\t".join(row)+ "\n")

					else:
						indef_shift.write("\t".join(row)+ "\n")


			else:
				indef.write(" ".join(row)+ "\n")

				if shift == ['NO']:
					indef_non_shift.write("\t".join(row)+ "\n")

				else:
					indef_shift.write("\t".join(row)+ "\n")


	U2_U12_like.close()
	indef.close()
	indef_shift.close()
	indef_non_shift.close()	
	
					

if __name__ == '__main__':
	main(sys.argv[1])	
