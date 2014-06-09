import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def main(sam):
	reader = csv.reader(open(sam), delimiter = '\t')
	
	inverse_flag = {"99":"147", "147":"99", "83":"163", "163":"83"  }
	
	for row in reader:
		if row[0] == "@SQ":
			print "\t".join(row)
		else:
			flag = row[1]
			rev_flag = inverse_flag[flag]
			seq = row[9]
			phred = row[10]
			rev_seq = str(Seq(seq).reverse_complement())
			rev_phred = row[10][::-1]			
			print row[0] +"\t" + rev_flag + "\t" + "\t".join(row[2:9]) + "\t" + rev_seq + "\t" + rev_phred + "\t" +  "\t".join(row[11:])

if __name__ == '__main__':
	main(sys.argv[1])
