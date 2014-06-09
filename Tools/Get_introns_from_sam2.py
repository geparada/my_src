import sys
import csv
import pysam


def main (SAM):
	
	""" Extrae los intrones de un SAM """
	samfile = pysam.Samfile( SAM, "r" )
	
	for row in samfile:
		start = row.pos + row.qstart
		
		print row.positions
	



if __name__ == '__main__':
	main(sys.argv[1])
	
