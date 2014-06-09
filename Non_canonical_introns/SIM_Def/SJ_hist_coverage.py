import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq

def main(SJ_check):
	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	reader2 = csv.reader(open(SJ_check), delimiter = ' ')	
		
	canonical_OK = defaultdict(int)
	canonical_wrong = defaultdict(int)
	non_canonical_OK = defaultdict(int)
	non_canonical_wrong = defaultdict(int)
	
	hist_canonical_OK = defaultdict(int)
	hist_canonical_wrong = defaultdict(int)
	hist_non_canonical_OK = defaultdict(int)
	hist_non_canonical_wrong = defaultdict(int)
	
	
	list_introns = []
	
	print >> sys.stderr, "Generando hush table con los intrones OK ..."

	for row in reader1:
		
		i = row[1]
		read_status = row[4]
		if read_status == 'OK':
			list_introns.append((i, 0))
	
	dict_introns = dict(list_introns)	
	
		
	print >> sys.stderr, "Contantdo numero de pruebas por intron ..."		
	
	introns_finded = []
	
	for row in reader2:
		i = row[1]
		read_status = row[4]
		dn = row[5]
		

		if read_status == 'OK' or read_status == 'Wrong':

			if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
				if dict_introns.has_key(i) == True:
					canonical_OK[i] += 1

				else:
					canonical_wrong[i] += 1

			else:
				if dict_introns.has_key(i) == True:
					non_canonical_OK[i] += 1

				else:
					non_canonical_wrong[i] += 1

			
			introns_finded.append(i)
			
	print >> sys.stderr, "Calculando estadisticas de numero de pruebas por intron ..."
			
	for i in set(introns_finded):
		hist_canonical_OK[canonical_OK[i]] += 1
		hist_canonical_wrong[canonical_wrong[i]] += 1
		hist_non_canonical_OK[non_canonical_OK[i]] += 1
		hist_non_canonical_wrong[non_canonical_wrong[i]] += 1
		
		
	print 'Coverage |OK Canonical | Wrong Canonical | %True discovery rate Canonical | OK Non-Canonical | Wrong Non-Canonical | %True discovery rate Non-Canonical |  '
		
	for n in range(1,1001):
		
		try:
			can_OK = hist_canonical_OK[n]
		except KeyError:
			can_OK = 0
		
		try:
			can_wrong = hist_canonical_wrong[n]
		except KeyError:
			can_wrong = 0
			
		try:	
			non_can_OK = hist_non_canonical_OK[n]
		except KeyError:
			non_can_OK = 0
		
		try:	
			non_can_wrong = hist_non_canonical_wrong[n]
		except KeyError:
			non_can_wrong = 0	

		print n, can_OK, can_wrong, percent(can_OK, can_OK+can_wrong), non_can_OK, non_can_wrong, percent(non_can_OK, non_can_OK+non_can_wrong)

			
def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0					
				


if __name__ == '__main__':
	main(sys.argv[1])  
