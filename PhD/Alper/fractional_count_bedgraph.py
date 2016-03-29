#!/nfs/users/nfs_g/gp7/miniconda/bin/python

import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict




def main(bam, chrominfo, target_chrom):

	bamfile = pysam.Samfile(bam, "rb")

        bedgraph_plus = open('bedgraph.plus.' + target_chrom, 'w')
        bedgraph_minus = open('bedgraph.minus.' + target_chrom, 'w')

        for row in  csv.reader(open(chrominfo), delimiter = '\t'):
                 chrom, size = row
                 size = int(size)


                 if chrom==target_chrom:



                        for n in range(size):

                                fractional_count_plus = float(0)
                                fractional_count_minus = float(0)

                                repeat = False

                                for read in bamfile.fetch(chrom, n, n+1):


                                        if read.is_reverse == False:

                                                NH = read.tags[0][1]
                                                fractional_count_plus += float(1)/ float(NH)


                                        else:

                                                NH = read.tags[0][1]
                                                fractional_count_minus += float(1)/ float(NH)


                                out_plus = "\t".join(map(str, [chrom, n, n+1, fractional_count_plus]))
                                out_minus = "\t".join(map(str, [chrom, n, n+1, fractional_count_minus]))


                                bedgraph_plus.write(out_plus + "\n")
                                bedgraph_minus.write(out_minus + "\n")


	bamfile.close()




if __name__ == '__main__':            
	 main(sys.argv[1], sys.argv[2], sys.argv[3])