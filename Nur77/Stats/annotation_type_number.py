import sys
import csv
from collections import defaultdict

Annotation_type_pie_chart = defaultdict(int)
Distance_to_TSS_hist = defaultdict(int)



def main(TF_FINAL_TABLE):

	reader = csv.reader(open(TF_FINAL_TABLE), delimiter = '\t')
	reader.next()

	for row in reader:


		PeakID, chrom, chromStart, chromEnd, strand, Peak_Score, Focus_Ratio__Region_Size, Annotation, Detailed_Annotation, Distance_to_TSS, Nearest_PromoterID, Entrez_ID, Nearest_Unigene, Nearest_Refseq, Nearest_Ensembl, Gene_Name, Gene_Alias, Gene_Description, Gene_Type, TF_ChromHMM, TF_DNAaseq, TF_enhancers, TF_motifs, TF_clusters = row
		
		TF_ChromHMM = TF_ChromHMM.split("-and-")
		TF_DNAaseq = TF_DNAaseq.split("-and-")
		TF_enhancers = TF_enhancers.split("-and-")
		TF_motifs = TF_motifs.split("-and-")
		TF_clusters = TF_clusters.split("-and-")

		Distance_to_TSS = int(Distance_to_TSS)

		Annotation_type_pie_chart[Annotation.split("(")[0]] += 1
		Distance_to_TSS_hist[Distance_to_TSS] += 1

		



	distance_range = 10000
	distance_hist_bar = 100



	for i in range(-distance_range/distance_hist_bar,distance_range/distance_hist_bar +1):
		
		int_ranges = []
		accomulated_freq_int = 0

		for n in range(i * distance_hist_bar -distance_hist_bar,i * distance_hist_bar):

			int_ranges.append(n)
			accomulated_freq_int += Distance_to_TSS_hist[n]

		print min(int_ranges), max(int_ranges), accomulated_freq_int





if __name__ == '__main__':
	main(sys.argv[1])