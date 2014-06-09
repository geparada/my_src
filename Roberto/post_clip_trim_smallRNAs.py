import sys
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

def main(fastq):
	
	for record in SeqIO.parse(fastq, "fastq"):
		seq = record.seq
		ID = record.id
		L = len(seq)
		kmer_len = 6
		Q = record.letter_annotations["phred_quality"]		
		
		if "AAAAAA" in record.seq and L >=100:
			
			i = 0
			start = i
			end = i + kmer_len 			
			seq_hex = str(seq[start:end])
			
			while i <= L:
				if seq_hex == "AAAAAA":
					break
				else:
					i += 1
					start = i
					end = i + kmer_len 					
					seq_hex = str(seq[start:end])										

			seq = seq[:i]
			Q = record.letter_annotations["phred_quality"][:i]
			
			
			if len(seq) >= 22:
					
				fastq_out = SeqRecord( seq, id = ID, description = "" )
				fastq_out.letter_annotations["phred_quality"] = Q
				
				print fastq_out.format("fastq"),
		

		else:
			fastq_out = SeqRecord( seq, id = ID, description = "" )
			fastq_out.letter_annotations["phred_quality"] = Q
				
			print fastq_out.format("fastq"),			

if __name__ == '__main__':
	main(sys.argv[1])  
