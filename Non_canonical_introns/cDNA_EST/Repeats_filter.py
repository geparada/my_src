import sys
import csv
from collections import defaultdict



def main(RepeatMasker, dust, introns):
	
	repeats = {}
	introns_ID = defaultdict(list)
	introns_keys = set([])
	
	reader1 = csv.reader(open(RepeatMasker), delimiter = '\t')
	reader2 = csv.reader(open(dust), delimiter = ' ')
	reader3 = csv.reader(open(introns), delimiter = ' ')
	
	for row in reader1:
		ID = row[9].split("|")[0]
		intron = row[9].split("|")[1]
		repeats[intron] = ID
	
	for row in reader2:
		ID = row[0][1:].split("|")[0]
		ID = row[0][1:].split("|")[1]
		repeats[intron] = ID
	
	for row in reader3:
		ID =  row[0]
		chr = row[1]
		istart = row[2]
		iend = row[3]
		strand = row[4]
		ilength = row[5]
		intron =row[6]
		dn = row[7]
		
		if repeats.has_key(intron)==False:
			introns_ID[intron].append(ID)
			introns_keys.add((chr, istart, iend, strand, ilength, intron, dn))

	
	for row in introns_keys:
		chr = row[0]
		istart = row[1]
		iend = row[2]
		strand = row[3]
		ilength = row[4]
		intron =row[5]
		dn = row[6]
		
		IDs = introns_ID[intron]
		coverage = len(IDs)
		
		print intron, coverage, chr, strand, istart, iend, ilength, dn, ",".join(IDs)


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])
