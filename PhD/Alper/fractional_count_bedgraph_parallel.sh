#!/bin/bash

chormInfo="/lustre/scratch108/compgen/team218/gp7/Genome/ce10/ce10_pEM975.chrom.sizes"

while getopts i:c: opt; do
	case $opt in
		i) Index=$OPTARG;;
		c) chrom=$OPTARG;;

	esac
done

name=$(echo $Index.$chrom)


echo $name


~/submit.job -n ${name}.process -s ~/my_src/PhD/Alper/fractional_count_bedgraph.py -q long -m 1000 $Index.sorted.bam $chormInfo $chrom 




