import sys
import csv

def main(file):

	for row in csv.reader(open(file), delimiter = '.'):
		print row[0]


if __name__ == '__main__':
	main(sys.argv[1])