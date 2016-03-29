import csv
import sys
import HTSeq
from collections import defaultdict 


def main(intronic_GTF, exonic_GTF, ENSEMBL_GTF):

	exon_out=open('WS220_gurdon_curated_exons_genes_pEM975_filtered.gtf', 'w+')
	intron_out=open('WS220_gurdon_curated_introns_genes_pEM975_filtered.gtf', 'w+')


	exons_Ensmbl = defaultdict(list)

	#It's necessary filter all the genes that won't have any introns alfter the filter

	gen_intron_count = defaultdict(int)
	gen_intron_count_ban = defaultdict(int)
	intron_ban = set([])


	for row in csv.reader(open(ENSEMBL_GTF), delimiter = '\t'):

		chrom, annotation, feature, start, end, score, strand, cero, ID = row


		start = int(start)
		end = int(end)

		chrom =  "chr" + chrom


		if feature=="exon":

			exon = HTSeq.GenomicInterval( chrom, start, end, "." )

			exons_Ensmbl[chrom].append(exon)

			#print exon

	for exon in exons_Ensmbl["chrX"]:

		intron = HTSeq.GenomicInterval( "chrX",	15174966, 15176772, ".")

		intron

	for row in csv.reader(open(intronic_GTF), delimiter = '\t'):

		chrom, annotation, feature, start, end, score, strand, cero, ID = row



		start = int(start)
		end = int(end)


		if feature=="exon":
			gen_intron_count[ID] += 1   
			intron = HTSeq.GenomicInterval( chrom, start, end, "." )
			intron_overlap = False

			

			for exon in exons_Ensmbl[chrom]:

				if intron.overlaps(exon):
					intron_overlap = True

			if intron_overlap:

				gen_intron_count_ban[ID] += 1
				intron_ban.add(intron)


	for row in csv.reader(open(intronic_GTF), delimiter = '\t'):

		chrom, annotation, feature, start, end, score, strand, cero, ID = row

		start = int(start)
		end = int(end)

		if feature =="gene":

			if gen_intron_count[ID] - gen_intron_count_ban[ID] != 0:

				intron_out.write("\t".join(row)+"\n")


		if feature=="exon":
			intron = HTSeq.GenomicInterval( chrom, start, end, "." )

			if (intron in intron_ban)==False:

				intron_out.write("\t".join(row)+"\n")


	for row in csv.reader(open(exonic_GTF), delimiter = '\t'):

		chrom, annotation, feature, start, end, score, strand, cero, ID = row

		if gen_intron_count[ID] - gen_intron_count_ban[ID] != 0:

			exon_out.write("\t".join(row)+"\n")




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])



#['chrIII', 'CE10_exonic', 'gene', '9226031', '9228033', '.', '-', '.', 'gene_id T23G5.1;']
#chrIII	CE10_exonic	gene	9226031	9228033	.	-	.	gene_id T23G5.1;

#chrI	ce6_ensGene	CDS	8378301	8378420	0.000000	-	0	gene_id "T19A6.1b"; transcript_id "T19A6.1b";