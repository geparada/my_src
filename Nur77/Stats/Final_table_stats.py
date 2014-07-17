import sys
import csv
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn2
import numpy as np

out_Annotation_type_pie_chart = open('Annotation_type_pie_chart', 'w')
out_Distance_to_TSS_hist = open('Distance_to_TSS_hist', 'w')
out_TF_ChromHMM_type_ocurrences = open('TF_ChromHMM_type_ocurrences', 'w')
out_intersections = open('peak_intersections', 'w')

out_DNAaseq_T = open('DNAaseq_T', 'w')
out_enhancers_T = open('enhancers_T', 'w')
out_motifs_T = open('motifs_T', 'w')
out_clusters_T = open('clusters_T', 'w')
out_ChromHMM_enhancers_T = open('ChromHMM_enhancers_T', 'w')
out_promoter_TSS_T = open('promoter_TSS_T', 'w')

Annotation_type_pie_chart = defaultdict(int)
Distance_to_TSS_hist = defaultdict(int)
TF_ChromHMM_type_ocurrences = defaultdict(int)

TOTAL = set([])

DNAaseq_T = set([])
enhancers_T = set([])
motifs_T = set([])
clusters_T = set([])
ChromHMM_enhancers_T = set([])
promoter_TSS_T = set([])


def main(TF_FINAL_TABLE):

	reader = csv.reader(open(TF_FINAL_TABLE), delimiter = '\t')
	reader.next()

	for row in reader:


		PeakID, chrom, chromStart, chromEnd, strand, Peak_Score, Focus_Ratio__Region_Size, Annotation, Detailed_Annotation, Distance_to_TSS, Nearest_PromoterID, Entrez_ID, Nearest_Unigene, Nearest_Refseq, Nearest_Ensembl, Gene_Name, Gene_Alias, Gene_Description, Gene_Type, TF_ChromHMM, TF_DNAaseq, TF_enhancers, TF_motifs, TF_clusters = row

		ID_Chr_Start_End = "_".join((chrom, chromStart, chromEnd))

		TOTAL.add(ID_Chr_Start_End)
		
		TF_ChromHMM = TF_ChromHMM.split("-and-")
		TF_DNAaseq = TF_DNAaseq.split("-and-")
		TF_enhancers = TF_enhancers.split("-and-")
		TF_motifs = TF_motifs.split("-and-")
		TF_clusters = TF_clusters.split("-and-")

		Distance_to_TSS = int(Distance_to_TSS)


		Annotation_type_pie_chart[Annotation.split("(")[0]] += 1
		Distance_to_TSS_hist[Distance_to_TSS] += 1

		if Annotation.split("(")[0]=="promoter-TSS ":
			promoter_TSS_T.add(ID_Chr_Start_End)

		
		for i in TF_ChromHMM:
			if i!="none":

				ID, ChromHMM_type = i.split("|")
				TF_ChromHMM_type_ocurrences[ChromHMM_type] += 1

				if ChromHMM_type.split("_")[-1]=="Enhancer":
					ChromHMM_enhancers_T.add(ID_Chr_Start_End)


		if TF_DNAaseq!=['none']:
			DNAaseq_T.add(ID_Chr_Start_End)

		if TF_enhancers!=['none']:
			enhancers_T.add(ID_Chr_Start_End)

		if TF_motifs!=['none']:
			motifs_T.add(ID_Chr_Start_End)

		for i in TF_clusters:

			if i!="none":
				cluster_ID, number_TFs, TFs = i.split("|")
				number_TFs = int(number_TFs)

				if number_TFs >= 2:
					clusters_T.add(ID_Chr_Start_End)


	### Peak annotation sites ####

	for k in Annotation_type_pie_chart.keys():
		out = [ k, Annotation_type_pie_chart[k] ]
		out_Annotation_type_pie_chart.write("\t".join(map(str, out)) + "\n")

	### TSS distance distribution ###

	distance_range = 10000
	distance_hist_bar = 100

	for i in range(-distance_range/distance_hist_bar +1,distance_range/distance_hist_bar +1):
		
		int_ranges = []
		accomulated_freq_int = 0

		for n in range(i * distance_hist_bar -distance_hist_bar,i * distance_hist_bar):

			int_ranges.append(n)
			accomulated_freq_int += Distance_to_TSS_hist[n]

		out = [ min(int_ranges), max(int_ranges), accomulated_freq_int ]

		out_Distance_to_TSS_hist.write("\t".join(map(str, out)) + "\n")


	for k in TF_ChromHMM_type_ocurrences.keys():
		out = [ k, TF_ChromHMM_type_ocurrences[k] ]
		out_TF_ChromHMM_type_ocurrences.write("\t".join(map(str, out)) + "\n")


	### Venn diagram ###

	#venn2([ChromHMM_enhancers_T,  clusters_T], ('ChromHMM_enhancers_T', 'clusters_T'))

	#plt.show()

	#venn3([ChromHMM_enhancers_T,  clusters_T, DNAaseq_T], ('ChromHMM_enhancers_T', 'clusters_T', 'DNAaseq_T'))

	#plt.show()

	#venn3([ChromHMM_enhancers_T, enhancers_T, clusters_T], ('ChromHMM_enhancers_T', 'enhancers_T', 'clusters_T'))

	#plt.show()

	#venn3([ChromHMM_enhancers_T, promoter_TSS_T, clusters_T], ('ChromHMM_enhancers_T', 'promoter_TSS_T', 'clusters_T'))

	#plt.show()	

	#venn2([promoter_TSS_T,  enhancers_T], ('promoter_TSS_T', 'enhancers_T'))

	#plt.show()

	#venn3([enhancers_T, promoter_TSS_T, ChromHMM_enhancers_T], ('enhancers_T', 'promoter_TSS_T', 'ChromHMM_enhancers_T'))

	#plt.show()

#	venn3([enhancers_T, promoter_TSS_T, motifs_T], ('enhancers_T', 'promoter_TSS_T', 'motifs_T'))
	
#	plt.show()

#	venn3([ChromHMM_enhancers_T, motifs_T, clusters_T], ('ChromHMM_enhancers_T', 'motifs_T', 'clusters_T'))

#	plt.show()	

	### Stacked bar graph ###



	intersections = []

	intersections.append('\t'.join(('ChromHMM_enhancers', str(len(ChromHMM_enhancers_T)), str(len(TOTAL) - len(ChromHMM_enhancers_T))  )))

	intersections.append('\t'.join(('FANTOM_enhancers', str(len(enhancers_T)), str(len(TOTAL) - len(enhancers_T))  )))

	intersections.append('\t'.join(('DNAaseq', str(len(DNAaseq_T)), str(len(TOTAL) - len(DNAaseq_T))  )))

	intersections.append('\t'.join(('promoter_TSS', str(len(promoter_TSS_T)), str(len(TOTAL) - len(promoter_TSS_T))  )))

	intersections.append('\t'.join(('clusters', str(len(clusters_T)), str(len(TOTAL) - len(clusters_T))  )))

	intersections.append('\t'.join(('motifs', str(len(motifs_T)), str(len(TOTAL) - len(motifs_T))  )))

	out_intersections.write('\n'.join(intersections))

	for i in ChromHMM_enhancers_T:
		out_ChromHMM_enhancers_T.write(i + '\n')

	for i in DNAaseq_T:
		out_DNAaseq_T.write(i + '\n')

	for i in enhancers_T:
		out_enhancers_T.write(i + '\n')

	for i in promoter_TSS_T:
		out_promoter_TSS_T.write(i + '\n')

	for i in clusters_T:
		out_clusters_T.write(i + '\n')

	for i in motifs_T:
		out_motifs_T.write(i + '\n')




	out_Annotation_type_pie_chart.close()
	out_Distance_to_TSS_hist.close()
	out_TF_ChromHMM_type_ocurrences.close()
	out_intersections.close()

	out_ChromHMM_enhancers_T.close()
	out_DNAaseq_T.close()
	out_ChromHMM_enhancers_T.close()
	out_promoter_TSS_T.close()
	out_clusters_T.close()
	out_motifs_T.close()



if __name__ == '__main__':
	main(sys.argv[1])