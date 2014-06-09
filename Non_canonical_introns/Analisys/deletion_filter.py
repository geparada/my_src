import sys
import csv

def main(deletions_table, TOTAL_final_table):
	
	csv.field_size_limit(sys.maxsize)
	
	reader1 = csv.reader(open(deletions_table), delimiter = ' ')
	reader2 = csv.reader(open(TOTAL_final_table), delimiter = ' ')
	
	deletions = set([])
	
	for row in reader1:
		deletion = row[0]
		deletions.add(deletion)

	
	for row in reader2:
		
		intron = row[0]

		
		if (intron in deletions) == False:
			print " ".join(row)

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])
