import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(ref_fasta,out_fasta):
	
	for ref,out in zip( SeqIO.parse(ref_fasta, "fasta"), SeqIO.parse(out_fasta, "fasta")):
		
		fasta_out = SeqRecord( out.seq, id = ref.id, description = "" )
		print fasta_out.format("fasta"),

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
