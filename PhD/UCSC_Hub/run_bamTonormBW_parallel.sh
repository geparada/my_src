#!/bin/bash


Genome="ce10"
fetchChromSizes $Genome > $Genome.chromsizes

for i in $(ls *.bam)


	do

		name=$(basename $i .bam)

		sed 's/name/$name/g' ~/my_src/PhD/UCSC_Hub/bamTonormBW.sh > bamTonormBW.$name.sh

		chmod +x bamTonormBW.$name.sh

		bash ~/submit.job -s ./bamTonormBW.$name.sh -m 50000 -q normal -n bw.$name.norm

	done
