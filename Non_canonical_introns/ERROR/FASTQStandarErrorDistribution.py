import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq

MAX_READ_LENGTH = 100
MAX_Q = 100

def ErrorCal(fastqfile, mutation_freq_index):
	reader = csv.reader(open(mutation_freq_index), delimiter=' ' )
	nt_P = []
	readscount = 0	
	for row in reader:
		nt = row[0]      #Nucleotide position number          
		P = row[2]       # Error probability
		nt_P.append((nt, P))

	index = dict(nt_P)
	

	f = open(fastqfile)

	scuare_mean_dif_sums = [0 for i in range(MAX_READ_LENGTH)]
	basecount = scuare_mean_dif_sums[:]
	logQs = [10 ** (float(-Q) / 10) for Q in range(MAX_Q + 1)]
	


	for record in SeqIO.parse(f, "fastq"):
		readscount += 1

		for bp, Q in enumerate(record.letter_annotations["phred_quality"]):
			scuare_mean_dif_sums[bp] += (logQs[Q] - float(index[str(bp)]))**2
			basecount[bp] += 1.0
				

	nextbasecount = basecount[1:] + [(0)]
	

	for bp, data in enumerate(zip(scuare_mean_dif_sums, basecount, nextbasecount )):
		scuare_mean_dif_sum, basecount, nextbasecount= data
		standar_error = (scuare_mean_dif_sum / readscount)**(0.5) / (readscount**(0.5))
		readlength = basecount - nextbasecount

		print bp, int(basecount), index[str(bp)], standar_error, int(readlength)


if __name__ == '__main__':
    ErrorCal(sys.argv[1],sys.argv[2])

