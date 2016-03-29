import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def PWM_to_dict(file):
	reader = csv.reader(open(file), delimiter = '\t')
	header = reader.next()
	header_dict = {}
	col = 0
	
	matrix = {}
	
	for name in header:
		header_dict[name] = col
		col += 1
	
	A_frec = []
	C_frec = []
	G_frec = []
	T_frec = []
	
	for row in reader:
		A = row[header_dict["A"]]
		C = row[header_dict["C"]]
		G = row[header_dict["G"]]
		T = row[header_dict["T"]]
		
		A_frec.append(float(A))
		C_frec.append(float(C))
		G_frec.append(float(G))
		T_frec.append(float(T))
			
	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	
	return matrix
def PWM_score(seq, PWM_dict):

	PWM_max_score = 0
	index = 0
	for N in range(len(PWM_dict["A"])):
		PWM_max_score += max(PWM_dict['A'][index], PWM_dict['C'][index], PWM_dict['T'][index], PWM_dict['G'][index])
		index += 1


	PWM_score = 0
	index = 0
	for N in seq:
		try:
			frec = PWM_dict[N][index]
			PWM_score += frec
			index += 1
		except KeyError:
			pass

	percentual_PWM_score = percent(PWM_score, PWM_max_score)

	return percentual_PWM_score

def main(motif_bed, PWM):


	PWM_dict = PWM_to_dict(PWM)

	for row in csv.reader(open(motif_bed), delimiter = '\t'):

		chrom, start, end, name, score, strand = row

		start = int(start)
		end = int(end)

		motif_seq = Genome[chrom][start: end]

		if strand == '-':
			motif_seq = motif_seq.reverse_complement()

		motif_seq = str(motif_seq).upper()

		print "\t".join(map(str, [name, motif_seq, PWM_score(motif_seq, PWM_dict)] )) 

		# new_score = (PWM_score(motif_seq, PWM_dict) - 80)*50

		# print "\t".join(map(str, [chrom, start, end, name, new_score, strand]))






if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3])
