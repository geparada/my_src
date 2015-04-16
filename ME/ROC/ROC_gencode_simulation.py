import sys
import csv
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np


def xfrange(start, stop, step):
    while start < stop:
        yield start
        start += step

def print_roc(neg, pos, decimals, step):

	start = round(min(neg + pos), decimals)
	end = round(max(neg + pos), decimals)

	for i in xfrange(start, end, step):
		print i, 1 - stats.percentileofscore(neg, i)/100, 1 - stats.percentileofscore(pos, i)/100




def main(gencode_scores, simulation_scores):

	gencode_U2 = []
	gencode_vertebrates = []
	gencode_primates = []

	simulation_U2 = []
	simulation_vertebrates = []
	simulation_primates = []


	for row in csv.reader(open(gencode_scores), delimiter = ' '):

		chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates  = row

		gencode_U2.append(float(U2_score))
		gencode_vertebrates.append(float(mean_conservation_vertebrates))
		gencode_primates.append(float(mean_conservation_primates))


	for row in csv.reader(open(simulation_scores), delimiter = ' '):

		chr, estart, eend, strand, U2_score, mean_conservation_vertebrates, mean_conservation_primates = row

		simulation_U2.append(float(U2_score))
		simulation_vertebrates.append(float(mean_conservation_vertebrates))
		simulation_primates.append(float(mean_conservation_primates))


	# min_U2 = round(min(gencode_U2 + simulation_U2), 1)
	# max_U2 = round(max(gencode_U2 + simulation_U2), 1)

	# for i in xfrange(min_U2, max_U2, 0.1):

	# 	print i, 1 - stats.percentileofscore(simulation_U2, i)/100, 1 - stats.percentileofscore(gencode_U2, i)/100


	#print_roc(simulation_U2, gencode_U2, 1, 0.1)
	#print_roc(simulation_vertebrates, gencode_vertebrates, 2, 0.01)
	print_roc(simulation_primates, gencode_primates, 2, 0.01)


	#plt.plot(U2_TPR, U2_FPR)
	#savefig('foo.png')
	#plt.clf()







if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])	