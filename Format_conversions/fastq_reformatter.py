import sys
from Bio import SeqIO
from Bio.Seq import Seq

def reformatter(fastq):
	
        f = open(fastq)

        for record in SeqIO.parse(f, "fastq"):
			print record.id
			print record.seq
			print record.letter_annotations["phred_quality"]

	f.close()


#		mutated = SeqRecord( Seq("".join(seq), generic_dna), id = record.id, description = "" )
#		mutated.letter_annotations["phred_quality"] = Q
#		print mutated.format("fastq")

if __name__ == "__main__":
    reformatter(sys.argv[1])
