#!/nfs/users/nfs_g/gp7/miniconda/bin/python

import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict




def main(bam1, bam2, RepeatMasker):

	bamfile1 = pysam.Samfile(bam1, "rb")
        bamfile2 = pysam.Samfile(bam2, "rb")

        
        #fractional_count = float(0)


        for row in  csv.reader(open(RepeatMasker), delimiter = '\t'):

                
                
                bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id= row


                normal_count1 = 0

                for read in bamfile1.fetch(genoName, int(genoStart), int(genoEnd)):

                        normal_count1 += 1

                normal_count2 = 0
                fractional_count2 = float(0)
                uniq_count2 = 0

                for read in bamfile2.fetch(genoName, int(genoStart), int(genoEnd)):

                        normal_count2 += 1


                        NH = read.tags[0][1]
                        fractional_count2 += float(1)/ float(NH)

                        if NH == 1:
                                uniq_count2 += 1



                out = map(str, [genoName, genoStart,  genoEnd,  repName, repClass, repFamily, normal_count1, normal_count2, fractional_count2, uniq_count2])

                print "\t".join(out)

	bamfile.close()




if __name__ == '__main__':            
	 main(sys.argv[1], sys.argv[2], sys.argv[3])