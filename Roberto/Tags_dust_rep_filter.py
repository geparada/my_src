import csv
import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


list_repeats_IDs = []
dict_repeats_IDs = {}

def repeats(RepeatMasker_psl, dust_fa):
	
	reader = csv.reader(open(RepeatMasker_psl), delimiter = '\t')
	
	for record in SeqIO.parse(dust_fa, "fasta"):
		
		if "N" in record.seq:
			list_repeats_IDs.append((str(record.id),record.seq))
	
	for row in reader:
		list_repeats_IDs.append((row[9],row[0]))
		
		



def main(tags):
	
	for record in SeqIO.parse(tags, "fasta"):
		
		sequence = str(record.seq)

		
		upperseq = SeqRecord( Seq(sequence.upper()), id = record.id, description = "" )
		#upperseq.letter_annotations["phred_quality"] = Q

		
		if dict_repeats_IDs.has_key(record.id)==False and len(sequence)==100:

			print upperseq.format("fasta"), 		

if __name__ == '__main__':
	repeats(sys.argv[1],sys.argv[2] )
	dict_repeats_IDs = dict(list_repeats_IDs)
	main(sys.argv[3])  
