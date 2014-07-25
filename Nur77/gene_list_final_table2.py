import sys
import csv
from collections import defaultdict

gene_list = set([])

def main(TF_FINAL_TABLE):

	reader = csv.reader(open(TF_FINAL_TABLE), delimiter = '\t')
	reader.next()

	for row in reader:

		PeakID, chrom, chromStart, chromEnd, strand, Peak_Score, Focus_Ratio__Region_Size, Annotation, Detailed_Annotation, Distance_to_TSS, Nearest_PromoterID, Entrez_ID, Nearest_Unigene, Nearest_Refseq, Nearest_Ensembl, Gene_Name, Gene_Alias, Gene_Description, Gene_Type, TF_ChromHMM, TF_DNAaseq, TF_enhancers, TF_motifs, TF_FANTOM_promoter_gene, TF_FANTOM_promoter_ID, TF_FANTOM_promoter_link, TF_clusters = row
		
		chromStart = int(chromStart)
		chromEnd = int(chromEnd)

		TF_ChromHMM = TF_ChromHMM.split("-and-")
		TF_DNAaseq = TF_DNAaseq.split("-and-")
		TF_enhancers = TF_enhancers.split("-and-")
		TF_motifs = TF_motifs.split("-and-")
		TF_clusters = TF_clusters.split("-and-")


		if "promoter-TSS" in Annotation:

			gene_list.add(Gene_Name)

		if TF_DNAaseq!=["none"]:

			gene_list.add(Gene_Name)

		if TF_enhancers!=["none"]:

			gene_list.add(Gene_Name)


		for i in TF_clusters:

			if i!="none":
				cluster_ID, number_TFs, TFs = i.split("|")
				number_TFs = int(number_TFs)

				if number_TFs >= 2:

					gene_list.add(Gene_Name)

		Distance_to_FATOM_TSS_genes = []

		TSS_genes = defaultdict(int)


		for gene, ID, link in zip(TF_FANTOM_promoter_gene.split("-and-"), TF_FANTOM_promoter_ID.split("-and-"),  TF_FANTOM_promoter_link.split(" ")):

			TSS_chrom = ""
			TSS_chromStart = ""
			TSS_chromEnd = ""
			TSS_strand = ""

			if gene != "unknown":
				TSS_genes[gene] += 1

			# if link != "none":
			# 	a, TSS_strand = link.split("http://fantom.gsc.riken.jp/5/sstar/FFCP_PHASE1:Hg19::")[1].split(",")
			# 	TSS_chrom, b = a.split(":")
			# 	TSS_chromStart, TSS_chromEnd = b.split("..")

			# 	TSS_chromStart = int(TSS_chromStart)
			# 	TSS_chromEnd = int(TSS_chromEnd)

			# 	Distance_to_FATOM_TSS = 0

			# 	if TSS_strand == "+":

			# 		Distance_to_FATOM_TSS = TSS_chromEnd - chromStart

			# 	if TSS_strand == "-":

			# 		Distance_to_FATOM_TSS = chromEnd - TSS_chromStart

			# 	if gene != "unknown":

			# 		Distance_to_FATOM_TSS_genes.append((Distance_to_FATOM_TSS, gene))

		True_Gene_Name = ""

		if TF_FANTOM_promoter_gene != "none":

			if Gene_Name in TSS_genes:
				True_Gene_Name = Gene_Name

			elif TSS_genes != {}:
				True_Gene_Name = max(TSS_genes.items() , key = lambda t: t[1])[0]

				print True_Gene_Name, Gene_Name, chrom, chromStart, chromEnd





		#print chrom, chromStart, chromEnd, Distance_to_FATOM_TSS_genes, max(TSS_genes.items() , key = lambda t: t[1]), Gene_Name



			#print gene, ID, link


	#for gene in gene_list:
		#print gene


if __name__ == '__main__':
	main(sys.argv[1])