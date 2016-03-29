for i in $(ls -d *ed?):

	do	
		cd $i

		echo $i

		rm *001.fastq
		rm STAR.worm.results/$i.Aligned.out.sam

		gzip ${i}_R1.fastq && rm ${i}_R1.fastq
		gzip ${i}_R2.fastq && rm ${i}_R2.fastq

		Pre-processing/${i}_R1.fastq && rm ${i}_R1.fastq
		Pre-processing/${i}_R2.fastq && rm ${i}_R2.fastq


		cd ..

	done
