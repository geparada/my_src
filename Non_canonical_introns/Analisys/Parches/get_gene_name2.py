import sys
import csv
from collections import defaultdict



def log_fold(fold):
	if fold == 0:
		return 0
	else:
		return math.log(fold**-1,2)			

def main(genecode_bed12, Final_table ):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')
	reader2 = csv.reader(open(genecode_bed12), delimiter = '\t')

	
	introns_exons_to_gene_name = {}
	gene_start_end = defaultdict(list)
	genes = set([])
	
		
	for row in reader2:


		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		end = int(row[2])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]
		gene = row[4]
		gene_start_end[chr+strand].append((gene, start, end, chr, strand))

		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			introns_exons_to_gene_name[chr+strand+str(istart)] = gene
			introns_exons_to_gene_name[chr+strand+str(iend)] = gene
			intron = row[0] + ":" +  str(istart) + row[5] + str(iend)


	c = 0

	for row in reader1:
		
		#row = row[1:]
		
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
		non_canonical_filter = row[26].split(",")
		
		TOTAL_variants = intron_retention_exon + skipped_exons_names + alt_introns + shift
		
		keys = set([])
		
		gene = "unknown_gene"

	

		for i in TOTAL_variants:
			if i!="NO":
				alt_chr = i.split("|")[0].split(":")[0]
				alt_start = i.split("|")[0].split(":")[1].split(strand)[0]
				alt_end = i.split("|")[0].split(":")[1].split(strand)[0]
				keys.add(alt_chr+strand+alt_start)
				keys.add(alt_chr+strand+alt_end)

		if introns_exons_to_gene_name.has_key(chr+strand+str(istart)):
			gene = introns_exons_to_gene_name[chr+strand+str(istart)]
		
		elif introns_exons_to_gene_name.has_key(chr+strand+str(iend)):
			gene = introns_exons_to_gene_name[chr+strand+str(iend)]
			
		
		else:
			for g in gene_start_end[chr+strand]:
				name = g[0]
				gen_start = g[1]
				gen_end = g[2]
				if (gen_start < istart and istart < gen_end) or (gen_start < iend and iend < gen_end):
					gene = name

			ID = ""
			
			
			
		if "unknown_gene" in gene:
			if ID != chr + str(istart)[:4] + strand:
				ID = chr + str(istart)[:4] + strand     #Como los datos estan ordenados previamente, esto servira para destiguir intrones que pertenesen a un mismo gen desconocido o a un distinto
				c += 1
				
			gene = "unknown_gene" + str(c)
		
		else:	
			genes.add(gene)
			
		print gene + "\t" + "\t".join(row)
		
#	print " ".join(genes)
	
#chr19	11,546,268	11,561,783
	 

#?       chr19:11,559,774+11559975 chr19   +       11559774        11559975        201     AACC    U12_ATAC        46.8051078907   0       12      0
 #      4       0       0       0       0       CCGCCTCTGCCCCTTCAAGCTTGTCTCGCAGAAACCCAACCGCTCCACCACCGA,GCAGAAACCCAACCGCTCCACCACCGTGCGCCTCCTGTGCGGGAAAGAGACCATGGTGACCAGCACCACAGAGCCC,GCAGAAACCCAACCGCTCCACCACCGTGCGCCTCCTGTGCGGGAAAGAGACCATGGTGACCAGCACCACAGAGCCC,CGAATACGTCTACCGCCTCTGCCCCTTCAAGCTTGTCTCGCAGAAACCCAACCGCTCCACCACCGTGCGCCTCCTG,CCAACGAATACGTCTACCGCCTCTGCCCCTTCAAGCTTGTCTCGCAGAAACCCAACCGCTCCAC  5       NO      NO      NO      NO
  #    NO      NO      NO      NO

			
			
if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])
