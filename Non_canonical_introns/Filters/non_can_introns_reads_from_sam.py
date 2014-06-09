import sys
import csv

def main(sam, spliced_reads, introns):
	
	reader1 = csv.reader(open(sam), delimiter = '\t')
	reader2 = csv.reader(open(spliced_reads), delimiter = ' ')
	reader3 = csv.reader(open(introns), delimiter = ' ')
	
	introns = []
	reads_list = []
	
	for row in reader3:
		intron = row[0]
		dn = row[7]
		introns.append((intron,dn))
		
	introns_dict = dict(introns)
	
	#print introns_dict
	
	for row in reader2:
		read = row[0]
		intron = row[6]
		#print intron
		if introns_dict.has_key(intron):
			reads_list.append((read,intron))
	
	reads_dict = dict(reads_list)
	
	#print reads_dict
	
	for row in reader1:
		if row[0] == "@SQ":
			print "\t".join(row)
		elif reads_dict.has_key(row[0]):
			print "\t".join(row)
	
	
if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])
