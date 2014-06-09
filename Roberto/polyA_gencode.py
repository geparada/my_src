import sys
import pysam
import csv

def main (bam, gencode_polyA):


	bamfile = pysam.Samfile(bam, "rb")
#	out = pysam.Samfile("out_plus.bam", "wb", template=bamfile)

	reader = csv.reader(open(gencode_polyA), delimiter = '\t')
	
	headers = reader.next()
	
	for row in reader:
		bin_, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames = row

		#if name2!="polyA_signal":
		if 0==0:

			for read in bamfile.fetch(chrom, int(txStart), int(txEnd)):
				qname, start, end, seq_len, flag, chr = read.qname, read.pos, read.aend, read.qlen, read.flag, bamfile.getrname(read.tid)
				print bam, qname, chr, start, end, name2			

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])	
