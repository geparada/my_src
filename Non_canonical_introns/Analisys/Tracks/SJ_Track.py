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

def main(Final_table):

	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)		
	
	for row in csv.reader(open(Final_table), delimiter = '\t'):  
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = float(row[8])
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		tissues_coverage = int(row[16])
		n_tissues = row[17]
		tissues = row[18]
		red = 0
		green = 0 
		blue = 0
		
		
		ID = dn[:2] + "-" + dn[2:]  + "[" + str(bodymap_coverage + gm12878_coverage + tissues_coverage) + "]"

		start = istart-8
		end = iend+8

		blockCount = "2"
		blockSizes = "8,8"
		blockStarts = "0" + "," + str(iend-start)
		
		start = istart-8
		end = iend+8
		
		if dn!="GTAG" and dn!="ATAC" and dn!="GCAG":
			
			if dn_type == "U2_GTAG":
				if dn_type_score >= score_U2_GTAG :
					green = int(100 + ((dn_type_score-score_U2_GTAG)*155)/(100-score_U2_GTAG))
				else:
					red = int(100 + ((score_U2_GTAG-dn_type_score)*155)/(score_U2_GTAG-20))

				
			elif dn_type == "U12_ATAC":
				if dn_type_score >= score_U12_ATAC:
					green = int(100 + ((dn_type_score-score_U2_GTAG)*155)/(100-score_U2_GTAG))
				else:
					red = int(100 + ((score_U12_ATAC-dn_type_score)*155)/(score_U12_ATAC-20))
			
				
			elif dn_type == "U12_GTAG":
				if dn_type_score >= score_U12_GTAG:
					green = int(100 + ((dn_type_score-score_U2_GTAG)*155)/(100-score_U2_GTAG))
				else:
					red = int(100 +  ((score_U12_GTAG-dn_type_score)*155)/(score_U12_GTAG-20))

	
		
		RGB = ",".join([str(red),str(green),str(blue)]) 
		
		BED = [chr, str(start), str(end ), ID, "0", strand, str(start), str(end), RGB, blockCount, blockSizes, blockStarts]
		print "\t".join(BED) 




if __name__ == '__main__':
	main(sys.argv[1])


