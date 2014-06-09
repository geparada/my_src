import sys
import csv
csv.field_size_limit(sys.maxsize)

def main(Final_table, dn_row):
	
	non_can = 0
	can = 0
	
	for row in csv.reader(open(Final_table), delimiter = ' '):
		dn = row[dn_row-1]
#		coverage = int(row[1])
		
#		if coverage>=3: 
	
		if dn=="GTAG" or dn=="GCAG" or dn=="ATAC":
			can += 1
		else:
			non_can += 1
		
	print "canonicos", can
	print "no canonicos", non_can


if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2]))
