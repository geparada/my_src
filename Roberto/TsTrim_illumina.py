import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import gzip

with_6Ts = open(sys.argv[1].strip(".gz") + ".6T", 'w')
without_6Ts = open(sys.argv[1].strip(".gz") + ".6T.noTs", 'w')

def main(fastq):
	
	
	for record in SeqIO.parse(fastq, "fastq"):

		seq = record.seq
		ID = record.id
		Q = record.letter_annotations["phred_quality"]
	
		if str(record.seq[:6]) =="TTTTTT" and set(record.seq)!=set(['T']):

			with_6Ts_fastq_out = SeqRecord( seq, id = ID, description = "" )
			with_6Ts_fastq_out.letter_annotations["phred_quality"] = Q
			
			with_6Ts.write(with_6Ts_fastq_out.format("fastq"),)

			new_seq = record.seq

			while new_seq[0]=="T":
				new_seq = new_seq[1:]
			
			if len(new_seq)>=16:
				
				without_6Ts_fastq_out = SeqRecord( new_seq, id = ID, description = "" )
				without_6Ts_fastq_out.letter_annotations["phred_quality"] = Q[-len(new_seq):]
				
				without_6Ts.write(without_6Ts_fastq_out.format("fastq"),)			


if __name__ == '__main__':
	main(sys.argv[1])  
