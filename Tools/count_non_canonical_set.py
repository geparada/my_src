import sys
import csv
import re
csv.field_size_limit(sys.maxsize)

def main(Final_table, dn_row):
		
	
	non_can = set([])
	can = set([])
	
	for row in csv.reader(open(Final_table), delimiter = ' '):
		intron = row[0]
#		coverage = int(row[1])
		dn = row[dn_row-1]
		dn = dn.upper()
		chr, start, end = re.findall(r"[\w']+", intron)		
		
		
#		if coverage>=3: 
		if int(end) - int(start)>=6:
		
			if dn=="GTAG" or dn=="GCAG" or dn=="ATAC":
				can.add(intron)
			else:
				non_can.add(intron)
		
	print "canonicos", len(can)
	print "no canonicos", len(non_can)


if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2]))
