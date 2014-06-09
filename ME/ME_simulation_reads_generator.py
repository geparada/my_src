import sys
import csv
from collections import defaultdict
from random import randint, sample


def main(simulation_seeds, Final_table_highfold):

	intron_coverages = []

	micro_exon_length_groups = defaultdict(list)




	for row in csv.reader(open(Final_table_highfold), delimiter = '\t'):

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
		tissues_coverage = int(tissues_coverage)
		n_tissues = int(n_tissues)
		intron_retention_exon = intron_retention_exon.split(",")
		skipped_exons_names = skipped_exons_names.split(",")
		alt_introns = alt_introns.split(",")
		alt_no_skipper_introns = alt_no_skipper_introns.split(",")
		alt_skipper_introns = alt_skipper_introns.split(",")
		alt_exon_variant_introns = alt_exon_variant_introns.split(",")
		shift = shift.split(",")
		non_canonical_shift = non_canonical_shift.split(",")

		coverage = bodymap_coverage + tissues_coverage

		if coverage >= 32:

			intron_coverages.append(coverage / 32)


	for row in csv.reader(open(simulation_seeds), delimiter = ' '):

		intron, micro_exon_tag, micro_exon_seq, micro_exon_start, micro_exon_end, intron3, intron5, micro_exon_iend, micro_exon_istart, new_intron3, new_intron5 = row

		micro_exon_length = len(micro_exon_seq)

		micro_exon_length_groups[micro_exon_length].append(row)



	for n in range(25):

		for row in sample(micro_exon_length_groups[n+1], 1000):

			intron, micro_exon_tag, micro_exon_seq, micro_exon_start, micro_exon_end, intron3, intron5, micro_exon_iend, micro_exon_istart, new_intron3, new_intron5 = row
			intron_coverage = sample(intron_coverages , 1)[0]

			#if intron_coverage > 32:
			#	intron_coverage = intron_coverage / 32

			#print intron, intron_coverage


			for r in range(intron_coverage):


				read_start = randint(0 ,len(micro_exon_tag)-100)
				read_end = read_start + 100

				read = micro_exon_tag[read_start:read_end]

				ID_info = [intron, micro_exon_seq, micro_exon_iend, micro_exon_istart, intron_coverage, r+1]
				ID_info = map(str, ID_info)

				print "@" + ("_").join(ID_info)
				print read
				print "+"
				print "H"*100







if __name__ == '__main__':
	main (sys.argv[1], sys.argv[2]) 		

