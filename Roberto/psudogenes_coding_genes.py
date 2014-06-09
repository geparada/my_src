import sys
import csv


def main(coding_genes, gencode):

	gencode_pseudogenes = set([])

	for row in csv.reader(open(gencode), delimiter = '\t'):

		gencode_pseudogenes.add(row[0])

	for row in csv.reader(open(coding_genes), delimiter = '\t'):

		gen = row[0]

		pseudogenes = []

		for i in gencode_pseudogenes:


			try:

				if i[:(len(gen) + 1)] == gen + "P":

						pseudogenes.append(i)

			except IndexError:
				pass


		n_pseudos = len(pseudogenes)

		print gen, n_pseudos, ",".join(pseudogenes)






if __name__ == '__main__':
		main (sys.argv[1], sys.argv[2])














