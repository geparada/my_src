import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(SJ_introns):
	
	reader1 = csv.reader(open(SJ_introns), delimiter = ' ')
	
	reads_to_print = set([])
	
	for row in reader1:
		ID = row[0]
		dn = row[7]
		seq = row[14]
		
		if dn!='GTAG' and dn!='GCAG' and dn!='ATAC' and dn!='CTAC' and dn!='GTAT' and dn!='CTGC':
			reads_to_print.add((ID, seq))
			
	for read in reads_to_print:
		ID = read[0]
		seq = read[1]
		print ">" + ID
		print seq
		
			
if __name__ == '__main__':
	
	main(sys.argv[1]) 
