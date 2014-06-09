import sys
import csv
import pysam


def main(bam):

	samfile = pysam.Samfile(bam, "rb")
#	output = pysam.Samfile("proccesed_multimapping.bam", "wb", template=samfile)
	
	for read in samfile.fetch():
		is_spliced = False
		matches = 0
		mismatches = 0
		deletions = 0
		insertions = 0

		for cigar_tuple in read.cigar:
			op = cigar_tuple[0]
			N = cigar_tuple[1]
			
			if op == 0:
				matches += N
			
			if op == 1:
				deletions += 1
			
			if op == 2:
				insertions += 1
			
			if op == 8:
				mismatches += N
		
		score = matches - mismatches - deletions - insertions
		
		if read.is_read1:
		
			print read.qname, 100*(float(score)/len(read.seq))
				
#		if is_spliced:
#			output.write(read)


#El score se calcula como matches - mismatches - numero qinsert - numero tinsert

	splicedreads.close()
	samfile.close()



if __name__ == '__main__':
	main(sys.argv[1])

