import sys
import pysam


def main (bam, chr, start, end):
	""" Imprime los reads que estan un una ventana genomica """

	bamfile = pysam.Samfile(bam, "rb")
	out_plus = pysam.Samfile("out_plus.bam", "wb", template=bamfile)
	out_minus = pysam.Samfile("out_minus.bam", "wb", template=bamfile)

	for alignedread in bamfile.fetch(chr, start, end):
		if alignedread.pos >= start and alignedread.aend <= end:
			strand = ""
			for t in alignedread.tags:
				if t[0] == "XS":
					strand = t[1]
			
			if strand == "+":
				out_plus.write(alignedread)

			elif strand == "-":
				out_minus.write(alignedread)							


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))	
