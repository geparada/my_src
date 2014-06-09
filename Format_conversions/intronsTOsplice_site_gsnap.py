import sys
import csv

def main(introns):
		reader1 = csv.reader(open(introns), delimiter = ' ')
		
		for row in reader1:
			ide = row[0] + '.' + row[6]
			chr = row[1]
			start = row[2]
			end = row[3]
			strand = row[4]
			length = row[5]
			if strand == '+':
				print '>' + ide, chr + ':' + start + '..' + str(int(start) + 1), 'donor', length 
				print '>' + ide, chr + ':' + end + '..' + str(int(end) + 1), 'acceptor', length
			else:
				print '>' + ide, chr + ':' + str(int(end) + 1) + '..' + end, 'donor', length
				print '>' + ide, chr + ':' + str(int(start) + 1) + '..' + start , 'acceptor', length 
			

if __name__ == '__main__':
	main(sys.argv[1])
