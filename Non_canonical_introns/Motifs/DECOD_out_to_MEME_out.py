import csv
import sys
from collections import defaultdict


def main (DECOD_out):
	""" Transforma el output de DECOD al un formato MEME """
	
	A = []
	C = []
	G = []
	T = []
	
	
	for row in csv.reader(open(DECOD_out), delimiter = ' '):
		try:
		
			if row[0] == "A":
				A.append(row[1:])
			elif row[0] == "C":
				C.append(row[1:])
			elif row[0] == "G":
				G.append(row[1:])
			elif row[0] == "T":
				T.append(row[1:])
		
		except IndexError:
			pass
	
	print "MEME version 4"
	print ""
	print "ALPHABET= ACGT"
	print ""
	print "strands: +"
	print ""
	print "Background letter frequencies"
	print "A 0.237 C 0.253 G 0.293 T 0.216"
	print "" 	
	
	n = 0
	
	for a, c, g, t in zip(A, C, G, T):
		
		n += 1
		
		for i in range(8):
			if i == 0:
				print "MOTIF", n
				print "letter-probability matrix: alength= 4 w= 8 nsites= 17 E= 4.1e-001"
				
			print a[i].strip('[').strip(']'), c[i].strip('[').strip(']'), g[i].strip('[').strip(']'), t[i].strip('[').strip(']')
		print ""
			



if __name__ == '__main__':
	main(sys.argv[1])
