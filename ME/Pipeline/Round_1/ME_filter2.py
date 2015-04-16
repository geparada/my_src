import sys
import csv
from collections import defaultdict
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict

Genome = {}

def uniform(seq):
	return str(seq).upper()


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def main(ME_centric, snp, non_canonical_introns):
	insertions = {}
	non_canonical_SJs = {}


	for row in csv.reader(open(snp), dialect=csv.excel_tab):
		_bin, chrom, chromStart, chromEnd, name, score, strand, refNCBI, refUCSC, observed, molType, _class, valid, avHet, avHetSE, func, locType, weight, exceptions, submitterCount, submitters, alleleFreqCount, alleles, alleleNs, alleleFreqs, bitfields = row

		insertion_freq = 0

		for a, f in zip(alleles.strip(",").split(","),  alleleFreqs.strip(",").split(",")):

			f = float(f)

			if a!="-":
				insertion_freq += f

		if refNCBI == "-":

			pos = " ".join([chrom, chromStart])

			info = " ".join([strand, observed, name, str(insertion_freq)])

			insertions[pos] = info

	for row in csv.reader(open(non_canonical_introns), dialect=csv.excel_tab):

		gene, intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_type_score, bodymap_coverage, gm12878_coverage, hg19_cDNA_coverage, hg19_EST_coverage,  mm9_cDNA_coverage, mm9_EST_coverage, genecode_coverage, tissues_coverage, n_tissues, tissues_names, intron_retention_exon, skipped_exons_names, alt_introns, alt_no_skipper_introns, alt_skipper_introns, alt_exon_variant_introns, shift, non_canonical_shift = row
		
		dn_matches = 0

		for nt, x in zip(dn, "GTAG"):
			if nt == x:
				dn_matches += 1

		if dn_matches>=3:

			non_canonical_SJs[intron] = float(dn_type_score)

	for row in csv.reader(open(ME_centric), delimiter = '\t'):

		############ SNP filter ###########

		sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, total_max_U2_scores, total_max_mean_conservations_vertebrates, total_max_mean_conservations_primates, P_ME, total_ME = row

		len_micro_exon_seq_found = int(len_micro_exon_seq_found)
		total_max_U2_scores = float(total_max_U2_scores)
		sum_total_coverage = int(sum_total_coverage)


		snp = False
		homopolymer = False
		burge = False


		for i in total_SJs.split(","):  #Analizando las SJ host

			chrom, chromStart, chromEnd = re.findall(r"[\w']+", i)

			chromStart = int(chromStart)
			chromEnd = int(chromEnd)

			for n in range(6):

				pos_start = " ".join([chrom, str(chromStart-n)])
				pos_end = " ".join([chrom, str(chromEnd+n)])

				if insertions.has_key(pos_start):

					strand, observed, name, insertion_freq = insertions[pos_start].split(" ")

					if ")" in observed:
						snp_at_SJ.add(True)

					else:
						for o in observed.split("/"):
							if len(o) == len_micro_exon_seq_found and o!="-":
								snp = True

				elif insertions.has_key(pos_end):

					strand, observed, name, insertion_freq = insertions[pos_end].split(" ")

					if ")" in observed:
						snp_at_SJ.add(True)

					else:
						for o in observed.split("/"):
							if len(o) == len_micro_exon_seq_found and o!="-":
								snp = True

			# if snp:
			# 	print row

			########## Homopolymer ###########

			strand = "+"

			if "-" in i:
				strand = "-"

			L = 50

			exon5 = uniform(Genome[chrom][chromStart-L:chromStart])
			intron5 = uniform(Genome[chrom][chromStart:chromStart+L])
			intron3 = uniform(Genome[chrom][chromEnd-L:chromEnd])
			exon3 = uniform(Genome[chrom][chromEnd:chromEnd+L])		
				
			if strand == "-":
				exon3 = uniform(Genome[chrom][chromStart-L:chromStart].reverse_complement())
				intron3 = uniform(Genome[chrom][chromStart:chromStart+L].reverse_complement())
				intron5 = uniform(Genome[chrom][chromEnd-L:chromEnd].reverse_complement())
				exon5 = uniform(Genome[chrom][chromEnd:chromEnd+L].reverse_complement())

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


			if poly_max >= 5 and nt_poly_max == "".join(set(nt_poly_max)):
			 	homopolymer = True



			############# Burge ############

			NC_intron_5 = chrom + ":" + str(chromStart + len_micro_exon_seq_found) + strand + str(chromEnd)

			NC_intron_3 = chrom + ":" + str(chromStart) + strand + str(chromEnd - len_micro_exon_seq_found)


			if intron5[:len_micro_exon_seq_found]==micro_exon_seq_found: 

				if intron5[len_micro_exon_seq_found:2+len_micro_exon_seq_found]=="GT":

					burge = True

				if NC_intron_5 in non_canonical_SJs:

					if non_canonical_SJs[NC_intron_5] > total_max_U2_scores:

						burge = True
				
			if intron3[-len_micro_exon_seq_found:]==micro_exon_seq_found:

				if intron3[-len_micro_exon_seq_found-2:-len_micro_exon_seq_found]=="AG":

					burge = True

				if NC_intron_3 in non_canonical_SJs:

					if non_canonical_SJs[NC_intron_3] > total_max_U2_scores:

						burge = True


		if snp==False and homopolymer==False and burge==False and sum_total_coverage>=3:

			print "\t".join(row)

		# if snp:
		# 	print "\t".join(row) + "\t" + "SNP"

		# if homopolymer:
		# 	print "\t".join(row) + "\t" + "homopolymer"

		# if burge:
		# 	print "\t".join(row) + "\t" + "Burge"

		# if sum_total_coverage < 3:
		# 	print "\t".join(row) + "\t" + "coverage"



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2], sys.argv[3], sys.argv[4])	


#python ~/my_src/ME/Pipeline/Round_1/ME_filter2.py ~/db/genome/hg19.fa TOTAL.sam.row_ME.filter1.ME_centric ~/db/Variation/snp138Common.fix ../../../Non_canonical_introns/TOTAL/non_canonical | sort -k 1 -n -r