from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
my_dna = Seq("GTAG,GCTG,ATAC", generic_dna)
print my_dna
print my_dna.complement()
print my_dna.reverse_complement()
