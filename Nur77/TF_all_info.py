import sys
import csv
from collections import defaultdict


def main(annotatePeaks, ChromHMM, DNAaseq, enhancers, motifs, TF_clusters, FANTOM_promoters ):

	ChromHMM_TF = defaultdict(list)
	DNAaseq_TF = defaultdict(list)
	enhancers_TF = defaultdict(list)
	motifs_TF = defaultdict(list)
	clusters_TF = defaultdict(list)
	FANTOM_promoters_TF = defaultdict(list) 


	for row in csv.reader(open(ChromHMM), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak, chrom2, chromStart2, chromEnd2, name2, score2, strand2 = row
		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))
		ID_feature = chrom2 + ":" + chromStart2 + "-" + chromEnd2 + "|" + name2
		ChromHMM_TF[ID_Chr_Start_End].append(ID_feature)


	for row in csv.reader(open(DNAaseq), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak, chrom2, chromStart2, chromEnd2, name2, score2, strand2 = row
		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))
		ID_feature = chrom2 + ":" + chromStart2 + "-" + chromEnd2
		DNAaseq_TF[ID_Chr_Start_End].append(ID_feature)

	for row in csv.reader(open(enhancers), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak, chrom2, chromStart2, chromEnd2, name2, score2, strand2, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts = row
		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))
		ID_feature = chrom2 + ":" + chromStart2 + "-" + chromEnd2
		enhancers_TF[ID_Chr_Start_End].append(ID_feature)

	for row in csv.reader(open(motifs), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak, chrom2, chromStart2, chromEnd2, name2, score2, strand2 = row
		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))
		ID_feature = chrom2 + ":" + chromStart2 + "-" + chromEnd2
		motifs_TF[ID_Chr_Start_End].append(ID_feature)

	for row in csv.reader(open(TF_clusters), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak, chrom2, chromStart2, chromEnd2, name2, score2, strand2 = row
		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))
		ID_feature = chrom2 + ":" + chromStart2 + "-" + chromEnd2 + "|" + score2  + "|" + name2 
		clusters_TF[ID_Chr_Start_End].append(ID_feature)


	for row in csv.reader(open(FANTOM_promoters), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak, chrom2, chromStart2, chromEnd2, name2, score2, strand2 = row
		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))

		chromStart2 = int(chromStart2)
		chromEnd2 = int(chromEnd2)

		up = -2000
		down = -500

		if strand2 == "+":

			chromStart2 = chromStart2 - up
			chromEnd2 = chromEnd2 + down

		if strand2 == "-":

			chromStart2 = chromStart2 - down
			chromEnd2 = chromEnd2 + up

		chromStart2 = str(chromStart2)
		chromEnd2 = str(chromEnd2)

		ID_feature =  (name2, chrom2, chromStart2, chromEnd2, strand2)

		FANTOM_promoters_TF[ID_Chr_Start_End].append(ID_feature)

	reader = csv.reader(open(annotatePeaks), delimiter = '\t')
	reader.next()

	print "\t".join(["PeakID", "chrom", "chromStart", "chromEnd", "strand", "Peak_Score", "Focus_Ratio/Region_Size", "Annotation", "Detailed_Annotation", "Distance_to_TSS", "Nearest_PromoterID", "Entrez_ID", "Nearest_Unigene", "Nearest_Refseq", "Nearest_Ensembl", "Gene_Name", "Gene_Alias", "Gene_Description", "Gene_Type", "ChromHMM", "TF_DNAaseq", "FANTOM_enhancer",  "TF_motif", "TF_FANTOM_promoter_ID", "TSS_gene", "TSS_Distance", "TSS_ID", "TSS_link", "Stronger_TSS_Distance", "Stronger_TSS_ID", "Stronger_TSS_link", "ENCODE_TF_clusters"])

	for row in reader:

		

		PeakID, chrom, chromStart, chromEnd, strand, Peak_Score, Focus_Ratio__Region_Size, Annotation, Detailed_Annotation, Distance_to_TSS, Nearest_PromoterID, Entrez_ID, Nearest_Unigene, Nearest_Refseq, Nearest_Ensembl, Gene_Name, Gene_Alias, Gene_Description, Gene_Type = row
		chromStart = str(int(chromStart)-1)

		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))

		TF_ChromHMM = ChromHMM_TF[ID_Chr_Start_End]
		TF_DNAaseq = DNAaseq_TF[ID_Chr_Start_End]
		TF_enhancers = enhancers_TF[ID_Chr_Start_End]
		TF_motifs = motifs_TF[ID_Chr_Start_End]
		TF_clusters = clusters_TF[ID_Chr_Start_End]
		TF_FANTOM_promoters = FANTOM_promoters_TF[ID_Chr_Start_End]

		TF_FANTOM_promoter_gene = []
		TF_FANTOM_promoter_ID = []
		TF_FANTOM_promoter_link = []

		if TF_ChromHMM == []:
			TF_ChromHMM = ["none"]
		else:
			TF_ChromHMM = ["-and-".join(TF_ChromHMM)]

		if TF_DNAaseq == []:
			TF_DNAaseq = ["none"]
		else:
			TF_DNAaseq = ["-and-".join(TF_DNAaseq)]

		if TF_enhancers == []:
			TF_enhancers = ["none"]
		else:
			TF_enhancers = ["-and-".join(TF_enhancers)]

		if TF_motifs == []:
			TF_motifs = ["none"]
		else:
			TF_motifs = ["-and-".join(TF_motifs)]

		if TF_clusters == []:
			TF_clusters = ["none"]
		else:
			TF_clusters = ["-and-".join(TF_clusters)]


		TSS_genes = defaultdict(int)
		Distance_to_FATOM_TSS_genes = defaultdict(list)

		TSS_gene = "unknown"

		TSS_Distance = "none"
		TSS_ID = "none"		
		TSS_link = "none"

		Stronger_TSS_Distance = "none"
		Stronger_TSS_ID = "none"		
		Stronger_TSS_link = "none"		

		if TF_FANTOM_promoters == []:
			TF_FANTOM_promoter_gene.append("none")
			TF_FANTOM_promoter_ID.append("none")
			TF_FANTOM_promoter_link.append("none")

		else:

			for i in TF_FANTOM_promoters:

				ID, TSS_chrom, TSS_chromStart, TSS_chromEnd, TSS_strand = i
				link = TSS_chrom + ":" + TSS_chromStart + ".." + TSS_chromEnd + "," + TSS_strand

				TSS_gene = "unknown"

				if "p@" != ID[:2]:

					TSS_gene = ID.split("@")[1]
					TSS_genes[TSS_gene] += 1



				link = "http://fantom.gsc.riken.jp/5/sstar/FFCP_PHASE1:Hg19::" + link

				TF_FANTOM_promoter_gene.append(TSS_gene)
				TF_FANTOM_promoter_ID.append(ID)
				TF_FANTOM_promoter_link.append(link)

				TSS_chromStart, TSS_chromEnd, chromStart, chromEnd = map(int, [TSS_chromStart, TSS_chromEnd, chromStart, chromEnd])



			 	Distance_to_FATOM_TSS = 0

			 	if TSS_strand == "+":

			 		Distance_to_FATOM_TSS = TSS_chromEnd - chromStart



			 	if TSS_strand == "-":

			 		Distance_to_FATOM_TSS = chromEnd - TSS_chromStart

			 	if TSS_gene != "unknown":

			 		Distance_to_FATOM_TSS_genes[TSS_gene].append((Distance_to_FATOM_TSS, ID, link))		


			max_TSS_gene_ocurrence = 0


			TSS_gene_candidates = []

			Positive_TSS_gene_candidates_Distances = []
			Negative_TSS_gene_candidates_Distances = []


			if TSS_genes != {}:

				max_TSS_gene_ocurrence = max(TSS_genes.items() , key = lambda t: t[1])[1]

			for k in TSS_genes.items():

				TSS_gene, TSS_gene_ocurrence = k

				if TSS_gene_ocurrence == max_TSS_gene_ocurrence:

					TSS_gene_candidates.append(TSS_gene)


			for TSS_gene in TSS_gene_candidates:

				for info in Distance_to_FATOM_TSS_genes[TSS_gene]:

					Distance, ID, link = info

					if Distance >= 0:

						Positive_TSS_gene_candidates_Distances.append((TSS_gene, Distance, ID, link))

					else:
						Negative_TSS_gene_candidates_Distances.append((TSS_gene, Distance, ID, link)) 

			
			if Positive_TSS_gene_candidates_Distances != []:

				TSS_gene, TSS_Distance, TSS_ID, TSS_link = min(Positive_TSS_gene_candidates_Distances , key = lambda t: t[1])

				TSS_gene, Stronger_TSS_Distance, Stronger_TSS_ID, Stronger_TSS_link = min(Positive_TSS_gene_candidates_Distances , key = lambda t: int(t[2].split("@")[0].strip("p")))

			elif Negative_TSS_gene_candidates_Distances:

				TSS_gene, TSS_Distance, TSS_ID, TSS_link = max(Negative_TSS_gene_candidates_Distances , key = lambda t: t[1])

				TSS_gene, Stronger_TSS_Distance, Stronger_TSS_ID, Stronger_TSS_link = min(Negative_TSS_gene_candidates_Distances , key = lambda t: int(t[2].split("@")[0].strip("p")))



			TSS_Distance, Stronger_TSS_Distance = map(str, [TSS_Distance, Stronger_TSS_Distance])


			TF_FANTOM_promoter_gene = ["-and-".join(TF_FANTOM_promoter_gene)]
			TF_FANTOM_promoter_ID = ["-and-".join(TF_FANTOM_promoter_ID)]
			TF_FANTOM_promoter_link = [" ".join(TF_FANTOM_promoter_link)]


		output = row + TF_ChromHMM + TF_DNAaseq + TF_enhancers + TF_motifs + TF_FANTOM_promoter_ID + [TSS_gene, TSS_Distance, TSS_ID, TSS_link, Stronger_TSS_Distance, Stronger_TSS_ID, Stronger_TSS_link] + TF_clusters  

		print "\t".join(output)





#python ~/my_src/Nur77/TF_all_info.py nur77.annotatePeaks K562/nur77.ChromHMM K562/nur77.DNase K562/nur77.enhancers K562/nur77.nur77_motif K562/nur77.TF_clusters




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])