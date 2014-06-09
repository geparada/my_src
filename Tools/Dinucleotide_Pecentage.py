import sys
import csv
from collections import defaultdict

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main(introns_final_table):
	csv.field_size_limit(1000000000)

	reader1 = csv.reader(open(introns_final_table), delimiter = ' ')
		
	dn_type = defaultdict(int)
	Total = 0
	dns = []
		
	for row in reader1:
		intron = row[0]
		dn = row[7]
		dn_type[dn] += 1
		Total += 1
		dns.append(dn)
		
	print "TOTAL =", Total 
	print "Dinucleotide_TYPE", "Number", "%"
		
	dn_frec = dn_type.items()
	dn_frec.sort(key=lambda x: x[1])
		
	for i in reversed(dn_frec):
		dn = i[0]
		frec = i[1]
			
		print dn, frec, percent(frec, Total)
						
if __name__ == '__main__':
	main(sys.argv[1])
