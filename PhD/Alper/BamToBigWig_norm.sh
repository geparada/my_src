#!/bin/bash


 for i in $(ls *.bam)

	 do 

	 	name=$(basename "$i" .bam)
	 	echo $name

	 	# echo "bamToBed ..." 
	 	# bamToBed -i $name.bam > $name.bed && echo "OK"

	 	# echo "Bed to BedGraph ..."  
	 	# bedItemOverlapCount $TEAM/Celegans_genome_ENSEMBL/CE10/ce10.fa $name.bed chromSize=$TEAM/Celegans_genome_ENSEMBL/CE10/ce10_pEM975.chrom.sizes | sort -k1,1 -k2,2n > $name.BedGraph && echo "OK"

		echo "normailized_bedgrap.py ..." 
		python ~/my_src/PhD/Alper/normalized_bedgraph.py $name.BedGraph > $name.norm.BedGraph && echo "OK"

		echo "sort ..." 
		awk '$1!="pEM975"' $name.norm.BedGraph| sort -k1,1 -k2,2n > $name.sort.norm.BedGraph && echo "OK"

		echo "bedGraphToBigWig ..."
		bedGraphToBigWig $name.sort.norm.BedGraph $TEAM/Celegans_genome_ENSEMBL/CE10/ce10_pEM975.chrom.sizes $name.norm.bw && echo "OK"


	 done

