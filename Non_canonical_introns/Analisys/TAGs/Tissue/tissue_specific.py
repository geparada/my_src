import sys
import csv

introns = set([])

def final_table_reader(Final_table):


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

		introns.add(intron)


def main (tissiue_ratios_table):
	
	reader = csv.reader(open(tissiue_ratios_table), delimiter = ' ')

	next(reader)

	results = []

	for row in reader:


		gene, intron, dn, alt_type, alt_donor,  alt_aceptor, alt_intron, alt_dn, adipose, adrenal, brain, breast, colon, heart, kidney, liver, lung, lymph_node, ovary, prostate, skeletal_muscle, testes, thyroid,  white_blood_cells = row

		tissues = [adipose, adrenal, brain, breast, colon, heart, kidney, liver, lung, lymph_node, ovary, prostate, skeletal_muscle, testes, thyroid,  white_blood_cells]
		tisues_names = ["adipose", "adrenal", "brain", "breast", "colon", "heart", "kidney", "liver", "lung", "lymph_node", "ovary", "prostate", "skeletal_muscle", "testes", "thyroid",  "white_blood_cells"]

		tissues_ratios = []

		total_non_canonical_coverage = {}
		total_canonical_coverage = {}

		recidual_factor = 0.1
		coverage_factor = 20

		#if dn!="GTAG" and dn!="ATAC" and dn!="GCAG":  #and gene == "GMFB":

		if intron in introns:


			for t, t_name in zip(tissues, tisues_names):

				non_canonical_coverage, canonical_coverage = map(float, t.split("/"))
				ratio = non_canonical_coverage / (canonical_coverage + recidual_factor)

				total_non_canonical_coverage[t_name] = non_canonical_coverage 
				total_canonical_coverage[t_name] = canonical_coverage

				if non_canonical_coverage + canonical_coverage >= coverage_factor:

					tissues_ratios.append((t_name, ratio))

			if len(tissues_ratios) >= 2:

				max_ratio_tissue, max_ratio = max(tissues_ratios, key=lambda x:x[1])

				del total_canonical_coverage[max_ratio_tissue]
				del total_non_canonical_coverage[max_ratio_tissue]

				sum_non_canonical =  sum(x[1] for x in total_non_canonical_coverage.items()) 
				sum_canonical = sum(x[1] for x in total_canonical_coverage.items())

				total_sum_ratio = (sum_non_canonical / ( sum_canonical + recidual_factor ) )  #* ((sum_non_canonical + sum_canonical ) / coverage_factor) 

				tissue_specificity_ratio = max_ratio / (total_sum_ratio + recidual_factor)

				#if sum_non_canonical >= coverage_factor:

				if sum_canonical >= coverage_factor:

					results.append((gene, intron, dn, alt_type, alt_donor,  alt_aceptor, alt_intron, alt_dn,  tissue_specificity_ratio, max_ratio_tissue, adrenal, brain, breast, colon, heart, kidney, liver, lung, lymph_node, ovary, prostate, skeletal_muscle, testes, thyroid,  white_blood_cells))


	results.sort(key=lambda x: x[8], reverse=True)


	for i in results:
		print " ".join(map(str, i))





if __name__ == '__main__':
	final_table_reader(sys.argv[2])
	main(sys.argv[1])