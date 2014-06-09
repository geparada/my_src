import sys
import csv

from Bio import SeqIO

def desconcatenador(text):

	reader1 = csv.reader(open(text), dialect='excel-tab' )
		
        for row in reader1:
		name = row[9].split(" ")
		for n in name:
			print  row[1] + '\t' + n


if __name__ == '__main__':
        desconcatenador(sys.argv[1])

