import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(introns,fastq):
	reader1 = csv.reader(open(introns), delimiter = ' ')
	
	list_introns = []
	
	for row in reader1:
		read = row[0]
		dn = row[7]
		if dn != 'GTAG' and dn != 'GCAG' and dn != 'ATAC':

			list_introns.append((read, 0))

	dict_intron = dict(list_introns)
	
	for record in SeqIO.parse(fastq, "fastq"):

		fastq_out = SeqRecord( record.seq, id = record.id, description = "" )
		fastq_out.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"]

		if dict_intron.has_key(str(record.id))==True:
			print fastq_out.format("fasta"),		
		
				
				


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])  
