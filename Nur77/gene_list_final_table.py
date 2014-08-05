import sys
import csv
from collections import defaultdict

gene_list = set([])

def main(TF_FINAL_TABLE):

	reader = csv.reader(open(TF_FINAL_TABLE), delimiter = '\t')
	reader.next()

	for row in reader:

		PeakID, chrom, chromStart, chromEnd, strand, Peak_Score, Focus_Ratio__Region_Size, Annotation, Detailed_Annotation, Distance_to_TSS, Nearest_PromoterID, Entrez_ID, Nearest_Unigene, Nearest_Refseq, Nearest_Ensembl, Gene_Name, Gene_Alias, Gene_Description, Gene_Type, TF_ChromHMM, TF_DNAaseq, TF_enhancers, TF_motifs, TF_FANTOM_promoter_ID, TSS_gene, TSS_Distance, TSS_ID, TSS_link, Stronger_TSS_Distance, Stronger_TSS_ID, Stronger_TSS_link, TF_clusters = row
	
		TF_ChromHMM = TF_ChromHMM.split("-and-")
		TF_DNAaseq = TF_DNAaseq.split("-and-")
		TF_enhancers = TF_enhancers.split("-and-")
		TF_motifs = TF_motifs.split("-and-")
		TF_clusters = TF_clusters.split("-and-")


#		if "promoter-TSS" in Annotation:
#			gene_list.add(Gene_Name)

		if TSS_gene != "unknown":
			gene_list.add(TSS_gene)

#		if TF_DNAaseq!=["none"]:
#			gene_list.add(Gene_Name)

#		if TF_enhancers!=["none"]:
#			gene_list.add(Gene_Name)

		if Annotation != "Intergenic":

			for i in TF_ChromHMM:

				if i!="none":

					ChromHMM_choords, ChromHMM_type = i.split("|")

					if ChromHMM_type=="4_Strong_Enhancer" or ChromHMM_type=="5_Strong_Enhancer":

						gene_list.add(Gene_Name)


			for i in TF_clusters:

				if i!="none":
					cluster_ID, number_TFs, TFs = i.split("|")
					number_TFs = int(number_TFs)

					if number_TFs >= 2:

						gene_list.add(Gene_Name)


	for gene in gene_list:
		print gene


if __name__ == '__main__':
	main(sys.argv[1])