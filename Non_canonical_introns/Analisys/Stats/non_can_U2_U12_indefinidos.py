import sys
import csv
from collections import defaultdict



def main(Final_table):
	reader1 = csv.reader(open(Final_table), delimiter = '\t')
	
	

	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)
	
	U2_GTAG = []
	U12_ATAC = []
	U12_GTAG = []
	indef = []

	for row in reader1:
		
#		print row[:17]
		
		gene, intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage, mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage = row[:17]

		#intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage, mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage = row[:16]


		
		istart = int(istart)
		iend = int(iend)
		ilength = int(ilength)
		dn_type_score = float(dn_type_score)
		bodymap_coverage = int(bodymap_coverage)
		gm12878_coverage = int(gm12878_coverage)
		hg19_cDNA_coverage = int(hg19_cDNA_coverage)
		hg19_EST_coverage  = int(hg19_EST_coverage)
		mm9_cDNA_coverage = int(mm9_cDNA_coverage)
		mm9_EST_coverage = int(mm9_EST_coverage)
		genecode_coverage = int(genecode_coverage)
		

		if dn!="GTAG" and dn!="ATAC" and dn!="GCAG":		
		
			if dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG and ilength >= 40:
				U2_GTAG.append(row)
				
			elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC and ilength >= 40:
				U12_ATAC.append(row)		
				
			elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG and ilength >= 40:
				U12_GTAG.append(row)
				
			else:
				indef.append(row)

			
#	print "dn_type", "number"
		
#	print "U2_GTAG"
#	for i in U2_GTAG:
#		print "\t".join(i)

#	print "U12_ATAC"
#	for i in U12_ATAC:
#		print "\t".join(i)		

#	print "U12_GTAG"
#	for i in U12_GTAG:
#		print "\t".join(i)		

#	print "indef"
	for i in  indef:
		print "\t".join(i)
	

				
	
	   	
			
					

if __name__ == '__main__':
	main(sys.argv[1])	
