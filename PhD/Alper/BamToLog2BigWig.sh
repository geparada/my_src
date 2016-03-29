#!/bin/bash

#samtools merge - ../Index12.T40.bam ../Index4.T40.bam ../Index8.T40.bam | samtools sort -o - WT_sorted > WT.bam &
#samtools merge - ../Index1.T40.bam ../Index5.T40.bam ../Index9.T40.bam | samtools sort -o - emb4a_sorted > emb4a.bam &
#samtools merge - ../Index10.T40.bam ../Index2.T40.bam ../Index6.T40.bam | samtools sort -o - emb4b_sorted > emb4b.bam &
#samtools merge - ../Index11.T40.bam ../Index3.T40.bam ../Index7.T40.bam | samtools sort -o - hrde1_sorted > hrde1.bam &


 # for i in $(ls *.bam)

 # do 
 # 	name=$(basename "$i" .bam)
 # 	echo $name
 # 	bamToBed -i $name.bam > $name.bed
 # 	bedItemOverlapCount $TEAM/Celegans_genome_ENSEMBL/CE10/ce10.fa $name.bed chromSize=$TEAM/Celegans_genome_ENSEMBL/CE10/ce10_pEM975.chrom.sizes | sort -k1,1 -k2,2n > $name.BedGraph

 # done


for i in $(echo emb4a emb4b hrde1)

	do 
		echo $i

		echo "Log2_bedgraph.py ..." 
		python ~/my_src/PhD/Alper/Log2_bedgraph.py WT.BedGraph $i.BedGraph $TEAM/Celegans_genome_ENSEMBL/CE10/ce10_pEM975.chrom.sizes > $i.log2.BedGraph && echo "OK"

		echo "Log2_bedgraph.py ..." 
		awk '$1!="pEM975"' $i.log2.BedGraph| sort -k1,1 -k2,2n > $i.sort.log2.BedGraph && echo "OK"

		echo "bedGraphToBigWig ..."
		bedGraphToBigWig $i.sort.log2.BedGraph $TEAM/Celegans_genome_ENSEMBL/CE10/ce10_pEM975.chrom.sizes $i.log2.bw && echo "OK"

	done


# python ~/my_src/PhD/Alper/Log2_bedgraph.py WT.BedGraph emb4a.BedGraph > emb4a.log2.BedGraph 
# python ~/my_src/PhD/Alper/Log2_bedgraph.py WT.BedGraph emb4b.BedGraph > emb4b.log2.BedGraph 
# python ~/my_src/PhD/Alper/Log2_bedgraph.py WT.BedGraph hrde1.BedGraph > hrde1.log2.BedGraph


# bedGraphToBigWig accepted_hits.plus.bedGraph ChromInfo.txt accepted_hits.plus.bw

#bamToBed -i accepted_hits.bam -split > accepted_hits.bed

#bedItemOverlapCount mm9 -chromSize=ChromInfo.txt | sort -k1,1 -k2,2n > accepted_hits.plus.bedGraph

#"emb4b-3", "hrde1-3", "WT-3", "emb4a-1", "emb4b-1", "hrde1-1", "WT-1", "emb4a-2", "emb4b-2", "hrde1-2", "WT-2", "emb4a-3"
#"./Index10.T40.bam" "./Index11.T40.bam" "./Index12.T40.bam" "./Index1.T40.bam" "./Index2.T40.bam" "./Index3.T40.bam" "./Index4.T40.bam" "./Index5.T40.bam" "./Index6.T40.bam" "./Index7.T40.bam" "./Index8.T40.bam" "./Index9.T40.bam"


#WT = 12, 4, 8
#emb4a = 1, 5, 9
#emb4b = 10, 2, 6
#hrde1 = 11, 3, 7