import sys
import csv

def main(final_table):
	reader1 = csv.reader(open(final_table), delimiter = ' ')
	
	for row in reader1:
		
		intron = row[0]
		coverage = row[1]
		chr = row[2]
		strand = row[3]
		start = row[4]
		end = row[5]
		length = row[6]
		dn =  row[7]
		dn_type = row[8]
		dn_type_score = row[9]		
		reads = row[10]
		BED = [chr, start, end, intron + "|" + dn_type + "|"  + dn_type_score + "|" + reads, coverage, strand, length, "0", dn]
		
		print "\t".join(BED)




if __name__ == '__main__':
	main(sys.argv[1]) 
