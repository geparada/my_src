#!/bin/bash

for i in $(ls Index* -d)

	do

		echo $i

		cd $i

			if [ ! -d "RSeQC_plots" ]; then mkdir RSeQC_plots; fi

				cd RSeQC_plots
		 			junction_saturation.py -i ../$i.Aligned.out.sam.C_SJ.sorted.bam -r $TEAM/Celegans_genome_ENSEMBL/CE10/ce10_Ensembl.bed12 -o $i
					infer_experiment.py -r $TEAM/Celegans_genome_ENSEMBL/CE10/ce10_Ensembl.bed12-i ../$i.Aligned.out.sam.C_SJ.sorted.bam > $i.infer

				cd ..

		cd ..

	done



#junction_saturation.py -i /lustre/scratch108/compgen/team218/ec12/Analysis/3_20150519/100worms/100worms.Aligned.out.sam.sorted.bam -r $TEAM/Celegans_genome_ENSEMBL/Caenorhabditis_elegans.WBcel235.81.bed12 -o 100worms
