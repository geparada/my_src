import sys
import csv
import pysam


def main(bam):

	samfile = pysam.Samfile(bam, "rb")
	splicedreads = pysam.Samfile("TOTAL.alignments.only-uniq.SJ.sort.bam", "wb", template=samfile)
	for read in samfile.fetch():
		is_spliced = False
		for cigar_tuple in read.cigar:
			if cigar_tuple[0]==3:
				is_spliced = True
				
		if is_spliced:
			splicedreads.write(read)

	splicedreads.close()
	samfile.close()



if __name__ == '__main__':
	main(sys.argv[1])

