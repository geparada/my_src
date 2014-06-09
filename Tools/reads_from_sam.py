import sys
import csv

def main(sam, spliced_reads):
	
	reader1 = csv.reader(open(sam), delimiter = '\t')
	reader2 = csv.reader(open(spliced_reads), delimiter = ' ')
	
	reads_list = []
	
	for row in reader2:
		reads_list.append((row[0],row[7]))
	
	reads_dict = dict(reads_list)
	
	for row in reader1:
		if row[0] == "@SQ":
			print "\t".join(row)
		elif reads_dict.has_key(row[0]):
			print "\t".join(row)
	
	
if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
