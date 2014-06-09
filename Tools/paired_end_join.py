import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(forward,reverse):
	
	F_list = []
	
	for F in SeqIO.parse(forward, "fastq"):	
		
		fastq_out_F = SeqRecord( F.seq, id = F.id, description = "" )
		fastq_out_F.letter_annotations["phred_quality"] = F.letter_annotations["phred_quality"]
		
		F_list.append( (F.id.strip("/1"), [fastq_out_F.format("fastq")]))
	
	F_dict = dict(F_list)
	
	for R in SeqIO.parse(reverse, "fastq"):			
			
		fastq_out_R = SeqRecord( R.seq, id = R.id, description = "" )
		fastq_out_R.letter_annotations["phred_quality"] = R.letter_annotations["phred_quality"]
		
		try:
			print F_dict[R.id.strip("/2")][0],
			print fastq_out_R.format("fastq"),
		
		except KeyError:
			pass
	
		

		
	

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])  
