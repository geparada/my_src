#!/bin/bash

Genome="ce10"




samtools sort -o  $name.sort.bam $name.bam
samtools index $name.sort.bam $name.sort.bam.bai

bamToBed -i $name.sort.bam -split > $name.bed

sort -k1,1 -k2,2n $name.bed | bedItemOverlapCount $Genome -chromSize=$Genome.chromsizes stdin  > $name.bedGraph

awk '$1!="sensor_piRNA_mjIs144"' $name.bedGraph > $name.bedGraph.filter

python ~/my_src/PhD/Alper/normalized_bedgraph.py $name.bedGraph.filter > $name.bedGraph.filter.norm

bedSort $name.bedGraph.filter.norm $name.bedGraph.filter.norm.sort 

bedGraphToBigWig $name.bedGraph.filter.norm.sort $Genome.chromsizes $name.bw





