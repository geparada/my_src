#!/nfs/users/nfs_g/gp7/miniconda/bin/python

import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict




def main(bam, RepeatMasker):

	bamfile = pysam.Samfile(bam, "rb")

        
        #fractional_count = float(0)


        for row in  csv.reader(open(RepeatMasker), delimiter = '\t'):

                normal_count = 0
                
                bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id= row

                for read in bamfile.fetch(genoName, int(genoStart), int(genoEnd)):


                        normal_count += 1


                        #NH = read.tags[0][1]
                        #fractional_count_minus += float(1)/ float(NH)

                out = map(str, [genoName, genoStart,  genoEnd,  repName, repClass, repFamily, normal_count])

                print "\t".join(out)

	bamfile.close()




if __name__ == '__main__':            
	 main(sys.argv[1], sys.argv[2])