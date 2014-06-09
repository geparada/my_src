import sys
import csv


def main(SNP, DVG):
	
	reader1 = csv.reader(open(SNP), delimiter = '\t')
	reader2 = csv.reader(open(DVG), delimiter = '\t')
	
	for row in reader1:
		
		chr = row[1]
		start = row[2]
		end = row[3]
		name = row[4]
		strand = row[6]
		type = row[11]	 # deletion	enum('unknown', 'single', 'in-del', 'het', 'microsatellite', 'named', 'mixed', 'mnp', 'insertion', 'deletion')	Class of variant (single, in-del, named, mixed, etc.)

		intron = chr + ":" + start + strand + end
		coverage = 0
		length = int(end) - int(start)
		dn = type
		reads = 0
		op_strand = ""
		if strand == "+":
			op_strand = "-"
		elif strand == "-":
			op_strand = "+"			

		if (type == "deletion" or type == "in-del") and length > 6:
			print intron, coverage, chr, strand, start, end, length, dn, reads
			print intron.replace(strand,op_strand) , coverage, chr, op_strand, start, end, length, dn, reads
			


	for row in reader2:
	

		chr = row[1]
		start = row[2]
		end = row[3]
		name = row[4]
		strand = row[6]
		type = row[11]   #	CopyNumber	varchar(255)	values	Type of variation
			
		intron = chr + ":" + start + strand + end
		coverage = 0
		length = int(end) - int(start)
		dn = type
		reads = 0
		if strand == "+":
			op_strand = "-"
		elif strand == "-":
			op_strand = "+"							

		if (type == "InDel") and length > 6:
			print intron, coverage, chr, strand, start, end, length, dn, reads
			print intron.replace(strand,op_strand) , coverage, chr, op_strand, start, end, length, dn, reads	
				

			
if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])
