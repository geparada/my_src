import sys
import csv


def main(FASTQ):
	
	F = open(FASTQ)
	
	for row in F:
		if len(row) >= 78:
			print row[:75]
		else:
			print row,
		
			
if __name__ == '__main__':
	main(sys.argv[1])  
