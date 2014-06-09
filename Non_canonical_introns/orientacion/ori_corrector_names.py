import sys
import csv
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

reori = []

def ori (file):
	reader = csv.reader(open(file), dialect='excel-tab' )
	for row in reader:
		DN = row[13].split(",")
		newrow = tuple(row)
		N = len(DN)
		if DN.count("GTAG") + DN.count("GCAG") - DN.count("CTAC") - DN.count("CTGC") < 0:
			newrow = row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], str(Seq(row[13], generic_dna).reverse_complement()) 
			
			if newrow[5].split(',')[0]=='+':

    				for c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 in zip(newrow[1].split(','), newrow[2].split(','), newrow[3].split(','), newrow[4].split(','), newrow[5].split(','), newrow[6].split(','), newrow[7].split(','), newrow[8].split(','), newrow[9].split(','), newrow[10].split(','), newrow[11].split(',')[::-1], newrow[12].split(',')[::-1], newrow[13].split(',')[::-1]):

					print c1, c2, c3, newrow[0], c4, '-', c1+c2+'-'+c3, c1+c2+'-'+c3+newrow[0], c8, c9, c10, c11, c12, c13
			else:
				for c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 in zip(newrow[1].split(','), newrow[2].split(','), newrow[3].split(','), newrow[4].split(','), newrow[5].split(','), newrow[6].split(','), newrow[7].split(','), newrow[8].split(','), newrow[9].split(','), newrow[10].split(','), newrow[11].split(',')[::-1], newrow[12].split(',')[::-1], newrow[13].split(',')[::-1]):

					print c1, c2, c3, newrow[0], c4, '+', c1+c2+'+'+c3, c1+c2+'+'+c3+newrow[0], c8, c9, c10, c11, c12, c13

 
		else:
			if newrow[5].split(',')[0]=='+':			

				for c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 in zip(newrow[1].split(','), newrow[2].split(','), newrow[3].split(','), newrow[4].split(','), newrow[5].split(','), newrow[6].split(','), newrow[7].split(','), newrow[8].split(','), newrow[9].split(','), newrow[10].split(','), newrow[11].split(','), newrow[12].split(','), newrow[13].split(',')):
		
					print c1, c2, c3, newrow[0], c4, c5, c6, c7, c8, c9, c10, c11, c12, c13

			else:
				for c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 in zip(newrow[1].split(','), newrow[2].split(','), newrow[3].split(','), newrow[4].split(','), newrow[5].split(','), newrow[6].split(','), newrow[7].split(','), newrow[8].split(','), newrow[9].split(','), newrow[10].split(','), newrow[11].split(','), newrow[12].split(','), newrow[13].split(',')):

					print c1, c2, c3, newrow[0], c4, c5, c6, c7, c8, c9, c10, c11, c12, c13




if __name__ == '__main__':
	ori(sys.argv[1])
