import sys
import csv

def main (SAM):
		reader1 = csv.reader(open(SAM), dialect='excel-tab' )
		for row in reader1:
			if row[0]!='@SQ':
				L = len(row[9])
				if L >= 50:
					print '\t'.join(row)
			else:
				print '\t'.join(row)
				



if __name__ == '__main__':
	main(sys.argv[1])  
