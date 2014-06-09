import sys
import csv

def main(final_table):
	reader1 = csv.reader(open(final_table), delimiter = ' ')
	
	for row in reader1:
		
		intron = row[0]
		chr = row[1]
		strand = row[2]
		start = row[3]
		end = row[4]
		length = row[5]
		dn =  row[6]
		
		hg19 = int(row[9])
		GM12878 = int(row[10])
		cDNA = int(row[11])
		EST = int(row[12])
		
		if hg19 >= 3 or GM12878 >= 3 or (cDNA + EST >= 3):	

			BED = [chr, start, end, "|".join(row), "0", strand, length, "0", dn]
			print "\t".join(BED)




if __name__ == '__main__':
	main(sys.argv[1]) 
