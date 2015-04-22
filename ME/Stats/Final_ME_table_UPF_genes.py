import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		transcript_SJ[chrfa.id.split("|")[1]] = chrfa.id.split("|")[0]

		
	print >> sys.stderr, "OK"

	f.close()

def main(Final_ME_TABLE, TAGs, gencode_bed12, UPF_genes):

	transcript_SJ = {}
	gene_transcript = {}
	shUPF_up_genes = set([])


	f = open(TAGs)
	for chrfa in SeqIO.parse(f, "fasta"):
		transcript_SJ[chrfa.id.split("|")[0]] = chrfa.id.split("|")[1]
	f.close()

	for row in csv.reader(open(gencode_bed12), delimiter = '\t'):
		
		csv.field_size_limit(1000000000)

		transcript = row[3]
		gene = row[4]

		gene_transcript[transcript] = gene


	for row in csv.reader(open(UPF_genes), delimiter = '\t'):

		if row[-1]=="yes" and ("-" in row[9])==True:
			shUPF_up_genes.add(row[2])

	for row in csv.reader(open(Final_ME_TABLE), delimiter = ' '):


		ME, total_SJs, U2_scores, mean_conservations_vertebrates, mean_conservations_primates, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, min_P_ME, score, GENCODE, Blencowe, Ponting, ME_cov_sum, SJ_cov_sum,  mixture_cov, adipose_cov, adrenal_cov, brain_cov, breast_cov, colon_cov, heart_cov, kidney_cov, liver_cov, lung_cov, lymph_node_cov, ovary_cov, prostate_cov, skeletal_muscle_cov, testes_cov, thyroid_cov, white_blood_cells_cov, HepG2_control_cov, HepG2_UPF2_cov, HELA_control_cov, HELA_UPF1_cov = row

		gene = gene_transcript[transcript_SJ[total_SJs.split(",")[0]]]

		if gene in shUPF_up_genes:

			print gene, ' '.join(row)



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


#python ~/my_src/ME/Stats/Final_ME_table_UPF_genes.py /media/HD3/Resultados/Micro_exons/Tags/Round2/Final/Final_ME_table /media/HD3/Resultados/Micro_exons/Tags/Round2/ME_canonical_SJ_tags.fa ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12