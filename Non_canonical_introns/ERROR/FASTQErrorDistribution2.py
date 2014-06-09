import sys
from Bio import SeqIO
from Bio.Seq import Seq

MAX_READ_LENGTH = 100
MAX_Q = 100

def ErrorCal(fastqfile):

	f = open(fastqfile)

	prob_sums = [0 for i in range(MAX_READ_LENGTH)]
	basecount = prob_sums[:]
	logQs = [10 ** (float(-Q) / 10) for Q in range(MAX_Q + 1)]
	


	for record in SeqIO.parse(f, "fastq"):

		for bp, Q in enumerate(record.letter_annotations["phred_quality"]):
			prob_sums[bp] += logQs[Q]
			basecount[bp] += 1.0	


	for bp, data in enumerate(zip(prob_sums, basecount)):
		prob_sum, basecount = data
		prob_mean = prob_sum / basecount
		print bp, int(basecount), prob_mean


if __name__ == '__main__':
    ErrorCal(sys.argv[1])

