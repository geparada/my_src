#!/bin/bash

Genome="ce10"
fetchChromSizes $Genome > $Genome.chromsizes

for i in $(ls *.bam)

	do

		name=$(basename $i .bam)


		samtools sort -o  $name.sort.bam $name.bam
		samtools index $name.sort.bam $name.sort.bam.bai

		bamToBed -i $name.sort.bam -split > $name.bed

		sort -k1,1 $name.bed | bedItemOverlapCount $Genome -chromSize=$Genome.chromsizes stdin  > $name.bedGraph

		awk '$1!="sensor_piRNA_mjIs144"' $name.bedGraph > $name.bedGraph.filter

		python ~/my_src/PhD/Alper/normalized_bedgraph.py $name.bedGraph.filter > $name.bedGraph.filter.norm

		bedGraphToBigWig $name.bedGraph.filter.norm $Genome.chromsizes $name.bw



	done

