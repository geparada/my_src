import sys
import csv

def main(file, c):

	for row in csv.reader(open(file), delimiter = '\t'):
		print row[c-1]


if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2]))