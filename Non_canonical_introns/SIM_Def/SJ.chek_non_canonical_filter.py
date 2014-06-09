import sys
import csv

def main(SJ_check):
	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	
	for row in reader1:
		i = row[1]
		status = row[4]
		dn = row[5]
		
		if dn != 'GTAG' and dn != 'GCAG' and dn != 'ATAC' and dn != "":
			print ' '.join(row)
		

if __name__ == '__main__':
	main (sys.argv[1])
