#!/bin/bash

for i in $(ls *.bam)


	do

		name=$(basename $i .bam)

		sed 's/name/$name/g' /Users/gp7/my_src/PhD/UCSC_Hub/bamTonormBW.sh > bamTonormBW.$name.sh

		chmod +x bamTonormBW.$name.sh

		bash /Users/gp7/submit.job -s ./bamTonormBW.$name.sh -m 50000 -q normal -n bw.$name.norm

	done
