import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable = []


def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = (str(chrfa.id), chrfa.seq)
		SeqTable.append(table)

	f.close()

     

def IntronExtractor(clip_trim):
         # row[0]  chr
         # row[1]  alignment start
         # row[2]  alignment end
         # row[3]  name (Refseq=intron0<>intron1=anchor0<>anchor1)
         # row[4]    
         # row[5]  strand
         # row[6]  aligment start
         # row[7]  aligment end
         # row[8]  
         # row[9] blocknum
         # row[10] blocksizes
         # row[11] qstarts  


	reader1 = csv.reader(open(clip_trim), delimiter = '\t')

	csv.field_size_limit(1000000000)	
	

	for row in reader1:
        	Genome = dict(SeqTable)

		#introns_set = set(introns_known)
		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))
		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]

		for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
			istart = start + q1 + b
			iend = start + q2
			ilen = iend - istart
			intron = row[0] +  str(istart) + row[5] + str(iend)

						
			dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]

			if row[5] == '-':
					dn = dn.reverse_complement()

			if ilen >= 40: #and str(dn).upper()!='GTAG' and str(dn).upper()!='GCAG' and str(dn).upper()!='ATAC' :						
				#print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (row[3], chr, istart, iend, strand, ilen, intron, str(dn).upper())
				print row[3], chr, istart, iend, strand, ilen, intron, str(dn).upper()
				

			#else:
			#	print row[3], chr, istart, iend, strand, ilen, intron, str(dn).upper()

				

			


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	IntronExtractor(sys.argv[2])

