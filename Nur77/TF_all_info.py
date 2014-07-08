import sys
import csv
from collections import defaultdict


def main(annotatePeaks, ChromHMM, DNAaseq, enhancers, motifs, TF_clusters  ):

	ChromHMM_TF = defaultdict(list)
	DNAaseq_TF = defaultdict(list)
	enhancers_TF = defaultdict(list)
	motifs_TF = defaultdict(list)
	clusters_TF = defaultdict(list)


	for row in csv.reader(open(ChromHMM), delimiter = '\t'):

		chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak, chrom2, chromStart2, chromEnd2, name2, score2, strand2 = row
		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))
		ID_feature = chrom2 + ":" + chromStart2 + "-" + chromEnd2 + "|" + name2
		ChromHMM_TF[ID_Chr_Start_End].append(ID_feature)

		#print ID_Chr_Start_End

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

	reader = csv.reader(open(annotatePeaks), delimiter = '\t')
	reader.next()

	print "\t".join(["PeakID", "chrom", "chromStart", "chromEnd", "strand", "Peak_Score", "Focus_Ratio/Region_Size", "Annotation", "Detailed_Annotation", "Distance_to_TSS", "Nearest_PromoterID", "Entrez_ID", "Nearest_Unigene", "Nearest_Refseq", "Nearest_Ensembl", "Gene_Name", "Gene_Alias", "Gene_Description", "Gene_Type", "ChromHMM", "TF_DNAaseq", "FANTOM_enhancer",  "TF_motif", "ENCODE_TF_clusters"])

	for row in reader:

		PeakID, chrom, chromStart, chromEnd, strand, Peak_Score, Focus_Ratio__Region_Size, Annotation, Detailed_Annotation, Distance_to_TSS, Nearest_PromoterID, Entrez_ID, Nearest_Unigene, Nearest_Refseq, Nearest_Ensembl, Gene_Name, Gene_Alias, Gene_Description, Gene_Type = row
		chromStart = str(int(chromStart)-1)

		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))

		TF_ChromHMM = ChromHMM_TF[ID_Chr_Start_End]
		TF_DNAaseq = DNAaseq_TF[ID_Chr_Start_End]
		TF_enhancers = enhancers_TF[ID_Chr_Start_End]
		TF_motifs = motifs_TF[ID_Chr_Start_End]
		TF_clusters = clusters_TF[ID_Chr_Start_End]

		if TF_ChromHMM == []:
			TF_ChromHMM = ["none"]

		if TF_DNAaseq == []:
			TF_DNAaseq = ["none"]

		if TF_enhancers == []:
			TF_enhancers = ["none"]

		if TF_motifs == []:
			TF_motifs = ["none"]

		if TF_clusters == []:
			TF_clusters = ["none"]


		output = row + TF_ChromHMM + TF_DNAaseq + TF_enhancers + TF_motifs + TF_clusters

		print "\t".join(output)





#python ~/my_src/Nur77/TF_all_info.py nur77.annotatePeaks K562/nur77.ChromHMM K562/nur77.DNase K562/nur77.enhancers K562/nur77.nur77_motif K562/nur77.TF_clusters




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])