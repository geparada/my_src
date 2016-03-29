import sys
import csv

def main(EISA_GTF):

	c = 0

	for row in csv.reader(open(EISA_GTF), delimiter = '\t'):

		c += 1

		chrom, source, feature_class, start, end, dot1, strand, dot2, IDs = row

		gene_name = IDs.split(" ")[1].strip(";")

		if feature_class == "gene":

			c = 0
			new_IDs = 'gene_id "' +  gene_name + '"'
			out = [chrom, source, "aggregate_gene", start, end, dot1, strand, dot2, new_IDs]

			print "\t".join(map(str, out))


		else:

			new_IDs = 'transcripts "' + gene_name + '"; exonic_part_number ' + '"' + str(c) + '"' + '; gene_id "' +  gene_name + '"'

			out = [chrom, source, "exonic_part", start, end, dot1, strand, dot2, new_IDs]

			print "\t".join(map(str, out))









if __name__ == '__main__':            
	 main(sys.argv[1])
