#!/bin/bash

Genome="ce10"




samtools sort NAME.bam NAME.sort
samtools index NAME.sort.bam NAME.sort.bam.bai

bamToBed -i NAME.sort.bam -split > NAME.bed

sort -k1,1 -k2,2n NAME.bed | bedItemOverlapCount $Genome -chromSize=$Genome.chromsizes stdin  > NAME.bedGraph

awk '$1!="sensor_piRNA_mjIs144"' NAME.bedGraph > NAME.bedGraph.filter

python ~/my_src/PhD/Alper/normalized_bedgraph.py NAME.bedGraph.filter > NAME.bedGraph.filter.norm

bedSort NAME.bedGraph.filter.norm NAME.bedGraph.filter.norm.sort 

bedGraphToBigWig NAME.bedGraph.filter.norm.sort $Genome.chromsizes NAME.bw





