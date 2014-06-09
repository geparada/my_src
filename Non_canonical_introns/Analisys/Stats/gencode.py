import sys
import csv
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main (gencode_introns, final_table):
	""" Calcula cuantos intrones estan en gencode """
	
	known_introns = set([])
	
	for row in csv.reader(open(gencode_introns), delimiter = ' '):
		
		name, chr, start, end, strand, lenght, intron, dn = row
		
		known_introns.add(intron)
	
	novel = 0
	known = 0
	
	conservated = 0
	
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
		
		if intron in known_introns:
			known += 1
		
		else:
			novel += 1

		if mm9_cDNA_coverage + mm9_EST_coverage >= 3:
			conservated += 1
	

	total =  novel + known

	print 'Nobel percentage'
	print novel, novel + known, percent(novel, total)

	print 'Conservation percentage'
	print conservated, total, percent(conservated, total)
	
		
	
	
	
if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])  
	
