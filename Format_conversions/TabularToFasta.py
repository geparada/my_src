import sys
import csv

def TTF(tabular):      
	reader1 = csv.reader(open(tabular), delimiter = ' ' )

	for row in reader1:
		print '>' + row[0] + '\n' + row[1]



if __name__ == '__main__':            
	 TTF(sys.argv[1])
