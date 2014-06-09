import sys
import csv
from Bio import SeqIO

def IDfinder(fasta,source,ID):
	
	f = open(fasta)
	for seq_record in SeqIO.parse(f, "fasta"):
		if source=="str":
			for seq_record in SeqIO.parse(f, "fasta"):
				if seq_record.id==ID:
					print ">" + seq_record.id + '\n' + seq_record.seq 
		if source=="list":
			csv_list = csv.reader(open(ID), delimiter=' ')
			for row in csv_list:
				try:
					IDs = row[4]
					if seq_record.id==IDs:
						print ">" + seq_record.id + '\n' + seq_record.seq		 			
				except IndexError:
					pass	
	f.close()


if __name__ == '__main__':                                      
	IDfinder(sys.argv[1],sys.argv[2],sys.argv[3])

