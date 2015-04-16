import sys
import csv
from collections import defaultdict 

def main(pre_processed, dust, repbase, genome_sam):

	read_SJ = defaultdict(set)
	black_list = set([])

	for row in csv.reader(open(dust), delimiter = '>'):

		black_list.add(row[1])

	for row in csv.reader(open(repbase), delimiter = '\t'):

		black_list.add(row[9])

	for row in csv.reader(open(genome_sam), delimiter = '\t'):

		if row[1]=="0" or row[1]=="16":

			black_list.add(row[0])

	for row in csv.reader(open(pre_processed), delimiter = '\t'):

		read, flag, tag, start, cigar, seq, qual = row

		SJ = tag.split("|")[0]
		read_SJ[read].add(SJ)


	for row in csv.reader(open(pre_processed), delimiter = '\t'):

		read, flag, tag, start, cigar, seq, qual = row

		if (read in black_list)==False and len(read_SJ[read])==1:
			print "\t".join(row)



		#print black_list


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])