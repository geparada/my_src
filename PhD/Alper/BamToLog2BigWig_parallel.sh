#!/bin/bash

Genome="$TEAM/Celegans_genome_ENSEMBL/CE10/ce10.fa"
Genome_chrom_sizes="/lustre/scratch108/compgen/team218/gp7/Genome/ce10/ce10_sensor_piRNA_mjIs144.chrom.sizes"
FILES="/nfs/users/nfs_t/tdd/lustre/projects/emb4/out/"


#samtools merge - $FILES/RNA_SX1316*.bam | samtools sort - WT &
#samtools merge - $FILES/RNA_SX2929*.bam | samtools sort - emb4a &
#samtools merge - $FILES/RNA_SX2930*.bam | samtools sort - emb4b &
#samtools merge - $FILES/RNA_SX2930*.bam | samtools sort - emb4 &
#samtools merge - $FILES/RNA_SX2000*.bam | samtools sort - hrde1 &


 for i in $(ls *.bam)

  do 
  	name=$(basename "$i" .bam)
  	echo $name
  	bamToBed -i $name.bam > $name.bed
  	bedItemOverlapCount $TEAM/Celegans_genome_ENSEMBL/CE10/ce10.fa $name.bed chromSize=$TEAM/Celegans_genome_ENSEMBL/CE10/ce10_pEM975.chrom.sizes | sort -k1,1 -k2,2n > $name.BedGraph

 done


for i in $(echo emb4 emb4a emb4b hrde1)

	do 
		echo $i

		echo "Log2_bedgraph.py ..."

		for chrom in $(awk '{print $1}' $Genome_chrom_sizes)

			do

				python ~/my_src/PhD/Alper/Log2_bedgraph_parallel.py WT.BedGraph $i.BedGraph $Genome_chrom_sizes $chrom > $i.log2.BedGraph.$chrom && echo $chrom

			done

		cat $i.log2.BedGraph.* > $i.log2.BedGraph && echo "OK"

		echo "Log2_bedgraph.py ..." 
		awk '$1!="pEM975"' $i.log2.BedGraph| sort -k1,1 -k2,2n > $i.sort.log2.BedGraph && echo "OK"

		echo "bedGraphToBigWig ..."
		bedGraphToBigWig $i.sort.log2.BedGraph $Genome_chrom_sizes $i.log2.bw && echo "OK"

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
