import sys
import csv
from collections import defaultdict
import numpy 

def main(TF_cluster, expNums_cell_file):

	expNums_cell_letter = defaultdict(set)

	for row in csv.reader(open(expNums_cell_file), delimiter = '\t'):

		expNums, letter, info1, info2, info3, info4, info5, info6 = row

		expNums_cell_letter[letter].add(expNums)


	for row in csv.reader(open(TF_cluster), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, expCount, expNums, expScores = row

		c = 0

		cell_expScores = []

		
		is_the_cell = False

		for i in expNums.split(','):

			if i in expNums_cell_letter['K']:  #Filtro para las K562

				is_the_cell = True
				cell_expScore = int(expScores.split(',')[c])
				cell_expScores.append(cell_expScore)

				c += 1

				
		score_cell = numpy.mean(cell_expScores)
		BED = [chrom, chromStart, chromEnd, name, score_cell, "."]

		

		if is_the_cell:
			print "\t".join(map(str, BED))



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])