import sys
import csv
from collections import defaultdict


def main (decamers_clusters, hexamers):
	""" Asigna cluster a hexameros de ESS """
	
	hexamer_cluster = set([])
	
	for a in csv.reader(open(hexamers), delimiter = '\t'):
		hexamer = a[0].upper()
		
		
		
		for b in csv.reader(open(decamers_clusters), delimiter = '\t'):
			decamer = b[0]
			cluster = b[1]

			if (hexamer in decamer) == True:
				hexamer_cluster.add((hexamer, cluster))
	
	for i in hexamer_cluster:
		print "\t".join(i)
		
	






if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])	
