
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC
from collections import defaultdict

Genome = {}

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0

def uniform(seq):
	return str(seq).upper()

def homopolymer (intron, chr, strand, istart, iend, Genome):
		L = 50

		exon5 = uniform(Genome[chr][istart-L:istart])
		intron5 = uniform(Genome[chr][istart:istart+L])
		intron3 = uniform(Genome[chr][iend-L:iend])
		exon3 = uniform(Genome[chr][iend:iend+L])		
			
		if strand == "-":
			exon3 = uniform(Genome[chr][istart-L:istart].reverse_complement())
			intron3 = uniform(Genome[chr][istart:istart+L].reverse_complement())
			intron5 = uniform(Genome[chr][iend-L:iend].reverse_complement())
			exon5 = uniform(Genome[chr][iend:iend+L].reverse_complement())

		poly5 = 0
		poly3 = 0
		poly_max = 0
		nt_poly_max = ''

		while exon5[-1] == exon5[-1 - poly5]:
			poly5 += 1
			if poly5 == L:
				break
		while exon3[0] == exon3[poly3]:
			poly3 += 1
			if poly3 == L:
				break

		if exon5[-1] == exon3[0]:
			poly_max = poly5 + poly3
			nt_poly_max = exon5[-1]
		else:
			poly_max, nt_poly_max = max((poly5, exon5[-1]), (poly3, exon3[0]), key=lambda x:x[0])

		return poly_max, nt_poly_max

def DR_counter(intron, ichr, strand, istart, iend, Genome):
	
	introns_finded_DR = [intron]
			
	L = 100         #Solo permite que se corra L pares de bases para buscar DR 		

	#Extrayendo regiones exonicas colindantes
			
	SJ5U = Genome[ichr][istart-L : istart].lower()
	SJ5D = Genome[ichr][istart : istart+L].lower()
	SJ3U = Genome[ichr][iend-L : iend].lower()
	SJ3D = Genome[ichr][iend : iend+L].lower()
	
	if strand == "-":
		SJ5U = Genome[ichr][iend : iend+L].lower().reverse_complement() 
		SJ5D = Genome[ichr][iend-L : iend].lower().reverse_complement()
		SJ3U = Genome[ichr][istart : istart+L].lower().reverse_complement()
		SJ3D = Genome[ichr][istart-L : istart].lower().reverse_complement()
	
	DRU = 0
	DRD = 0
								
	#Contando directos repetidos y generando intrones no consenso alternativos
		
	try:
		while SJ5U[L-1-DRU]==SJ3U[L-1-DRU]:
			DRU += 1
			if strand == "+":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart-DRU) + strand + str(iend-DRU)]
			elif strand == "-":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart+DRU) + strand + str(iend+DRU)]
			if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
				break
	except IndexError:
		pass 
	try:
		while SJ5D[DRD]==SJ3D[DRD]:
			DRD += 1
			if strand == "+":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart+DRD) + strand + str(iend+DRD)]
			elif strand == "-":
				introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart-DRD) + strand + str(iend-DRD)]
			if SJ5D[DRD]!=SJ3D[DRD]:
				break
	except IndexError:
		pass
	
	return DRU + DRD


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[str(chrfa.id)] = chrfa.seq

	f.close()  

	print >> sys.stderr, "OK"

def main(gencode_introns, Final_table):

	gencode_ID_introns = defaultdict(list)
	score_U2_GTAG = float(63.2891027914)
	score_U12_ATAC = float(60.9280810964)
	score_U12_GTAG = float(61.4553595446)

	header = ["Gene name", "Splice Junction ID", "Chromosome", "Strand", "Start", "End", "Intron length", "Dinucleotides", "Type-like", "score", "Mixture of 16 tissues reads", "GM12878 reads", "Human cDNAs", "Human ESTs", "Total 16 tissue samples reads", "Number of tissues", "Tissues",  "%GC", "Direct Repeats", "Homopolymer length", "Homopolymer nucleotide", "GENCODE v.17 ID"]

	print "\t".join(header)

	for row in csv.reader(open(gencode_introns), delimiter = ' '):
		
		name, chr, start, end, strand, lenght, intron, dn = row

		gencode_ID_introns[intron].append(name)



	for row in csv.reader(open(Final_table), delimiter = '\t'):
		
		gene, intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, tissues_coverage, n_tissues, tissues_names, intron_retention_exon, skipped_exons_names, alt_introns, alt_no_skipper_introns, alt_skipper_introns, alt_exon_variant_introns, shift, non_canonical_shift = row
		
		istart = int(istart)
		iend = int(iend)
		ilength = int(ilength)
		dn_type_score = float(dn_type_score)
		bodymap_coverage = int(bodymap_coverage)
		gm12878_coverage = int(gm12878_coverage)
		hg19_cDNA_coverage = int(hg19_cDNA_coverage)
		hg19_EST_coverage  = int(hg19_EST_coverage)
		mm9_cDNA_coverage = int(mm9_cDNA_coverage)
		mm9_EST_coverage = int(mm9_EST_coverage)
		genecode_coverage = int(genecode_coverage)
		n_tissues = int(n_tissues)
		intron_retention_exon = intron_retention_exon.split(",")
		skipped_exons_names = skipped_exons_names.split(",")
		alt_introns = alt_introns.split(",")
		alt_no_skipper_introns = alt_no_skipper_introns.split(",")
		alt_skipper_introns = alt_skipper_introns.split(",")
		alt_exon_variant_introns = alt_exon_variant_introns.split(",")
		shift = shift.split(",")
		non_canonical_shift = non_canonical_shift.split(",")

		intron_seq = str(Genome[chr][istart:iend]).upper()
		DR = DR_counter(intron, chr, strand, int(istart), int(iend), Genome)
		poly_max, nt_poly_max = homopolymer(intron, chr, strand, int(istart), int(iend), Genome)

		if poly_max == 1:
			nt_poly_max = "None"


		GC = intron_seq.count("G") + intron_seq.count("C")

		GC_percent = percent(GC, len(intron_seq))

		DR = DR_counter(intron, chr, strand, int(istart), int(iend), Genome)

		genecode_ID = "Unannotated"

		if intron in gencode_ID_introns:
			genecode_ID = ",".join(gencode_ID_introns[intron])



		if dn!="GTAG" and dn!="ATAC" and dn!="GCAG":

			min_ilength = 85
			high_score = 70

			rescue = alt_introns != ['NO'] or (mm9_EST_coverage + mm9_cDNA_coverage) >= 3

			if ilength >= min_ilength:		
		
				if dn_type == "U2_GTAG" and dn_type_score >= high_score:
					pass


				elif dn_type == "U2_GTAG" and dn_type_score >= score_U2_GTAG and rescue:
					pass
					
				elif dn_type == "U12_ATAC" and dn_type_score >= high_score:
					pass

				elif dn_type == "U12_ATAC" and dn_type_score >= score_U12_ATAC and rescue:
					pass

					
				elif dn_type == "U12_GTAG" and dn_type_score >= high_score:
					pass

				elif dn_type == "U12_GTAG" and dn_type_score >= score_U12_GTAG and rescue:
					pass
					
				else:
					dn_type = "Non-U2/U12"


			else:
				dn_type = "Non-U2/U12"


		output = [gene, intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage, tissues_coverage, n_tissues, tissues_names, GC_percent, DR, poly_max, nt_poly_max, genecode_ID ]

		print "\t".join(map(str, output))


		#print GC(intron_seq)



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3])

#python ~/my_src/Non_canonical_introns/Analisys/paper_final_table.py ~/db/genome/hg19.fa ~/db/transcriptome/hg19/Gene_models/gencode/v17/gencode.v17.annotation.introns non_canonical.highfold.indef