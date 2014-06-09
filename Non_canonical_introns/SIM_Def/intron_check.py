import sys
import csv
from collections import defaultdict

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main(row_introns, SJ_check):
	
	intron_OK_list = []
	Canonical = defaultdict(int)
	Non_canonical = defaultdict(int)
	
	Canonical_OK = []
	Non_canonical_OK = []
	
	Canonical_OK_and_no_found = []
	Non_canonical_OK_and_no_found = []
	
	reader1 = csv.reader(open(row_introns), delimiter = ' ')	
	reader2 = csv.reader(open(SJ_check), delimiter = ' ')
	
	for row in reader2:
		intron = row[1]
		status = row[4]
		dn = row[5]
		
		if status == 'NO_SJ_found' or status == 'OK':
			if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
				Canonical_OK_and_no_found.append(intron)
			else:
				Non_canonical_OK_and_no_found.append(intron)
		
			if status == 'OK':
				intron_OK_list.append((intron, status))
				if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
					Canonical_OK.append(intron)
				else:	
					Non_canonical_OK.append(intron)
						
	intron_OK_dict = dict(intron_OK_list)		
		
	Total_Canonical = len(set(Canonical_OK_and_no_found))
	Total_non_canonical = len(set(Non_canonical_OK_and_no_found))	
		
	for row in reader1:
		intron = row[0]
		dn = row[7]
		
		try:
			status = intron_OK_dict[intron]
			
		except KeyError:
			status = "Wrong"
	
		if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
			Canonical[status] += 1
		
		else:
			Non_canonical[status] += 1
	
	
	print 'Intron Type | Total | OK | Wrong | %True discovery rate | %Sensitivity'
	 	
	print 'Canonical', Total_Canonical, Canonical['OK'], Canonical['Wrong'], percent(Canonical['OK'], Canonical['OK'] + Canonical['Wrong']), percent(Canonical['OK'], Total_Canonical)
	print 'Non_canonical', Total_non_canonical, Non_canonical['OK'], Non_canonical['Wrong'], percent(Non_canonical['OK'], Non_canonical['OK'] + Non_canonical['Wrong']), percent(Non_canonical['OK'], Total_non_canonical)
			
		
		
	

if __name__ == '__main__':
	
	main(sys.argv[1],sys.argv[2]) 
