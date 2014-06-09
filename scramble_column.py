import sys
import csv
from random import shuffle


def shuffle_word(word):
	word = list(word)
	shuffle(word)
	return ''.join(word)

def shuffle_column(file):
	reader = csv.reader(open(file), dialect='excel-tab' )
	for row in reader:
		scramble_seq = shuffle_word(row[13])
		print row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], scramble_seq



if __name__ == '__main__':                                      #Permite tomar argumentos mediante el terminal
    shuffle_column(sys.argv[1])

