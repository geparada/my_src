import csv
import sys


def main(gene_featurecounts, exon_featurecounts):

	accepted_genes = set([])


	for row in csv.reader(open(gene_featurecounts), delimiter = '\t'):

		if row[0][0]!="#":

			geneID = row[0]
			accepted_genes.add(geneID)



	for row in csv.reader(open(exon_featurecounts), delimiter = '\t'):


		if row[0][0]=="#":

			print "\t".join(row)

		elif row[0] =="sensor_piRNA_mjIs144":
			print "\t".join(row)

		elif row[0] in accepted_genes:

			print "\t".join(row)










if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])