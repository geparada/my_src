import sys
import csv


def main(sam):
	reader = csv.reader(open(sam), delimiter = '\t')
	
	Positive_strand = set([])
	Negative_strand = set([])
	
	for row in reader:
		flag = row[1]
		
		if row[0] != "@SQ" and "N" in row[5]:
		
			if "XS:A:+" in row:
				Positive_strand.add(flag)
				
			elif "XS:A:-" in row:
				Negative_strand.add(flag)
			
	print "strand +", Positive_strand
	print "strand -", Negative_strand
	

if __name__ == '__main__':
	main(sys.argv[1])
		
		
		
