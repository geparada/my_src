import sys
import csv
import pysam
from decimal import *
getcontext().prec = 20
csv.field_size_limit(1000000000)

#M 0 alignment match (can be a sequence match or mismatch)
#I 1 insertion to the reference
#D 2 deletion from the reference
#N 3 skipped region from the reference
#S 4 soft clipping (clipped sequences present in SEQ)
#H 5 hard clipping (clipped sequences NOT present in SEQ)
#P 6 padding (silent deletion from padded reference)
#= 7 sequence match
#X 8 sequence mismatch

final_table_introns = set([])

def track_gen(sub_final_table, RGB):	
	
	
	for row in csv.reader(open(sub_final_table), delimiter = ' '):
		  
		intron = row[0]
		coverage = int(row[1])
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		ilength = row[6]
		dn = row[7]
		dn_type = row[8]
		dn_type_score = float(row[9])
		IDs = row[10].split(",")
		
		if intron in final_table_introns:
			
			for ID in IDs:
			
				start = istart-8
				end = iend+8

				blockCount = "2"
				blockSizes = "8,8"
				blockStarts = "0" + "," + str(iend-start)
				
				start = istart-8
				end = iend+8
				
				BED = [chr, str(start), str(end ), ID, "0", strand, str(start), str(end), RGB, blockCount, blockSizes, blockStarts]
				print "\t".join(BED) 


def main (mrna, est, Final_table):

	for row in csv.reader(open(Final_table), delimiter = '\t'):  
		intron = row[0]
		final_table_introns.add(intron)
	
	track_gen(mrna, "205,155,29")  #goldenrod3
	track_gen(est, "0,104,139")   #deepskyblue4		
	

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2], sys.argv[3])


