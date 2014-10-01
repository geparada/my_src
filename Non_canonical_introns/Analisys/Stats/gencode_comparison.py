import sys
import csv
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord



def main (old_gencode, new_gencode, final_table):
	""" Calcula cuantos intrones estan en gencode """
	
	old_gencode_introns = set([])
	new_gencode_introns = set([])
	our_introns = set([])	
	
	for row in csv.reader(open(old_gencode), delimiter = ' '):
		
		name, chr, start, end, strand, lenght, intron, dn = row
		
		old_gencode_introns.add(intron)

	for row in csv.reader(open(new_gencode), delimiter = ' '):
		
		name, chr, start, end, strand, lenght, intron, dn = row
		
		new_gencode_introns.add(intron)

	
	for row  in csv.reader(open(final_table), delimiter = '\t'):
		
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

		our_introns.add(intron)


	print (our_introns & old_gencode_introns) - new_gencode_introns
	print (our_introns & new_gencode_introns) - old_gencode_introns





		

	
if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])  
	
