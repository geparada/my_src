import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict


Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq.split(" ")[0]
		
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


def main(GTF, U2_GTAG_5_file, U2_GTAG_3_file, GFF_DEXSEQ):

	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)
		
	U2_GTAG_5_win = range(13)
	U2_GTAG_5_max_score = 0
	index = 0
	for N in U2_GTAG_5_win:
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])
		index += 1
	
	U2_GTAG_3_win = range(17)
	U2_GTAG_3_max_score = 0
	index = 0
	for N in U2_GTAG_3_win:
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])
		index += 1

	transcripts_exons = defaultdict(list)

	for row in csv.reader(filter(lambda row: row[0]!='#', open(GTF)), delimiter = '\t'):

		chrom, source, feature, start, end, score, strand, frame, attribute =  row

		
		if feature=='exon':

			transcript = attribute.split("; ")[1].split(" ")[1]
			transcripts_exons[transcript].append(row)



	donor_acceptor = {}
	acceptor = {}

	for i in transcripts_exons.items():

		transcript, rows = i

		if len(rows)>=3:

			for row in rows[1:-1]:

				chrom, source, feature, start, end, score, strand, frame, attribute =  row

				start = int(start)
				end = int(end)

				estart = start - 1
				eend = end

				exon_5 = Genome[chrom][(estart-14):(estart+3)]
				exon_3 = Genome[chrom][(eend-3):(eend+10)]


				if strand=="-":

					exon_5 = Genome[chrom][(eend-3):(eend+14)]
					exon_3 = Genome[chrom][(estart-10):(estart+3)]

					exon_5 = exon_5.reverse_complement()
					exon_3 = exon_3.reverse_complement()

				exon_5 = str(exon_5).upper()
				exon_3 = str(exon_3).upper()

				U2_GTAG_5_score = 0
				index = 0 
				for N in exon_3:
					try:
						frec = U2_GTAG_5[N][index]
						U2_GTAG_5_score += frec
						index += 1
					except KeyError:
						pass

				U2_GTAG_3_score = 0
				index = 0 
				for N in exon_5:
					try:
						frec = U2_GTAG_3[N][index]
						U2_GTAG_3_score += frec
						index += 1					 
					except KeyError:
						pass

				U2_GTAG_total_score = percent( (U2_GTAG_5_score + U2_GTAG_3_score), (U2_GTAG_5_max_score + U2_GTAG_3_max_score) )

				U2_GTAG_5_score = percent( U2_GTAG_5_score, U2_GTAG_5_max_score )
				U2_GTAG_3_score = percent( U2_GTAG_3_score, U2_GTAG_3_max_score )
				chrom = chrom.strip("chr")

				if strand == "+":

					donor_acceptor["_".join(map(str, [chrom, start, strand]))] = "\t".join([str(U2_GTAG_3_score), "acceptor"])
					donor_acceptor["_".join(map(str, [chrom, start-1, strand]))] = "\t".join([str(U2_GTAG_3_score), "acceptor"])

					donor_acceptor["_".join(map(str, [chrom, end, strand]))] = "\t".join([str(U2_GTAG_5_score), "donor"])
					donor_acceptor["_".join(map(str, [chrom, end+1, strand]))] = "\t".join([str(U2_GTAG_5_score), "donor"])


				if strand == "-":

					donor_acceptor["_".join(map(str, [chrom, end, strand]))] = "\t".join([str(U2_GTAG_3_score), "acceptor"])
					donor_acceptor["_".join(map(str, [chrom, end+1, strand]))] = "\t".join([str(U2_GTAG_3_score), "acceptor"])

					donor_acceptor["_".join(map(str, [chrom, start, strand]))] = "\t".join([str(U2_GTAG_5_score), "donor"])
					donor_acceptor["_".join(map(str, [chrom, start-1, strand]))] = "\t".join([str(U2_GTAG_5_score), "donor"])



	for row in csv.reader(open(GFF_DEXSEQ), delimiter = '\t'):


		chrom, gff_file, feature, start, end, dot1, strand, dot2, IDs = row

		start_score = []
		end_score = []

		if feature == "exonic_part":

				gene = IDs.split(" ")[-1].strip('"')

					

				if strand == "+":

					try:

						start_score =  donor_acceptor["_".join([chrom, start, strand])]

					except KeyError:
						pass

					try:

						end_score = donor_acceptor["_".join([chrom, end, strand])]

					except KeyError:
						pass



				if strand == "-":

					try: 

						start_score = donor_acceptor["_".join([chrom, end, strand])]

					except KeyError:
						pass 

					try:

						end_score = donor_acceptor["_".join([chrom, start, strand])]

					except KeyError:
						pass 



		if start_score!=[] and end_score!=[]:

			print "\t".join(row) + "\t" + start_score + "\t" + end_score



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]) # , sys.argv[2], sys.argv[2])


# 15 82109521 82109538 76.6534370656 []

#python ~/my_src/PhD/Teichmann/U2_score_DEXSeq.py GRCm38.p4.genome.fa gencode.vM9.chr_patch_hapl_scaff.annotation.gtf U2.GT_AG.donor.pwm U2.GT_AG.acceptor.pwm gencode.vM9.chr_patch_hapl_scaff.annotation.DEXseq.gff.sed
