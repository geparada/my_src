import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0		


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

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

def main(sim_fastq,  U2_GTAG_5_file, U2_GTAG_3_file, phylop_vertebrates, phylop_primates):

	MEs = set([])
	wWigIO.open(phylop_vertebrates)
	wWigIO.open(phylop_primates)

	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])
	
	TOTAL_U2_max_score = U2_GTAG_5_max_score + U2_GTAG_3_max_score

	for row in csv.reader(open(sim_fastq), delimiter = '\t'):


		if row[0][0]=="@":

			SJ, ME_seq, estart, eend, total_coverage, n = row[0].split("_")

			chr = SJ.split(":")[0][1:]
			strand = "+"

			if "-" in SJ:
				strand = "-"

			estart = int(estart)
			eend = int(eend)

			MEs.add((chr, strand, estart, eend))


	for m in MEs:

		chr, strand, estart, eend = m

		estart, eend = sorted([estart, eend])

		E5 = str(Genome[chr][estart-14:estart+3]).upper()
		E3 = str(Genome[chr][eend-3:eend+10]).upper()


		if strand == "-":

			E5 = str(Genome[chr][eend-3:eend+14].reverse_complement()).upper()
			E3 = str(Genome[chr][estart-10:estart+3].reverse_complement()).upper()


		U2_score = 0

		i = 0

		for N in E5:
			U2_score += U2_GTAG_3[N][i]
			i += 1

		i = 0

		for N in E3:
			U2_score += U2_GTAG_5[N][i]
			i += 1

		U2_score = percent(U2_score, TOTAL_U2_max_score)

		conservation_vertebrates = wWigIO.getIntervals(phylop_vertebrates, chr, estart-2, eend+2)
		conservation_primates = wWigIO.getIntervals(phylop_primates, chr, estart-2, eend+2)

		mean_conservation_vertebrates = 0
		mean_conservation_primates = 0

		for i in conservation_vertebrates:

			mean_conservation_vertebrates += i[2]

		try:

			mean_conservation_vertebrates = mean_conservation_vertebrates/len(conservation_vertebrates)

		except ZeroDivisionError:
			pass

		
		for i in conservation_primates:

			mean_conservation_primates += i[2]

		try:

			mean_conservation_primates = mean_conservation_primates/len(conservation_primates)

		except ZeroDivisionError:
			pass

		print E5, E3, chr, strand, estart, eend, U2_score, mean_conservation_vertebrates, mean_conservation_primates












if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main (sys.argv[2]) 		