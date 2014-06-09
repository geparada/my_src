import sys 
from Bio import SeqIO
from Bio.Seq import Seq




def Genomictabulator(fastqfile):
	
	f = open(fastqfile)

	Ptotal = []

	
	for record in SeqIO.parse(f, "fastq"):
		Pseq = []		

		for Q in record.letter_annotations["phred_quality"]:
			P = 10**(float(-Q)/10)
			Pseq = Pseq + [P]

		Ptotal = Ptotal + [Pseq]
	

	for l in range(100):
		Pn = []
		for P in Ptotal:
						
			try:	
				Pn = Pn + [P[l]]
			except IndexError:
				pass

		print l, len(Pn), float(sum(Pn))/float(len(Pn))				
	

		



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])   
