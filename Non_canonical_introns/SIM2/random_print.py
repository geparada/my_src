import sys
import csv
import random

def ran_printer(bigfile, P):
	
	reader = csv.reader(open(bigfile), dialect='excel-tab' )

	for row in reader:
		R = random.random()
		if R<=P:
			print " ".join(row)

if __name__ == '__main__':
	ran_printer(sys.argv[1], float(sys.argv[2]) )
