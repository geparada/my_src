
files=$(echo wgEncodeCshlShortRnaSeqA549CytosolShorttotalCiptapRawDataRep3.fastq wgEncodeCshlShortRnaSeqA549CytosolShorttotalTapRawDataRep3.fastq wgEncodeCshlShortRnaSeqA549NucleusShorttotalCiptapRawDataRep3.fastq wgEncodeCshlShortRnaSeqA549NucleusShorttotalTapRawDataRep3.fastq)



function preprocessing () {

	name=${1%.fastq}
	echo $name

	echo "Pre-processing $name..." $(date) >> $name.temp.log

	cat $name.fastq | echo $name.fastq $((`wc -l`/4)) >> $name.temp.log
	fastx_clipper -a AAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -l 22 -i $name.fastq -o $name.fastq.clip || fastx_clipper -a AAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -l 22 -i $name.fastq -o $name.fastq.clip -Q 33

	cat $name.fastq.clip | echo $name.fastq.clip $((`wc -l`/4)) >> $name.temp.log
	python ~/my_src/Roberto/post_clip_trim_smallRNAs.py $name.fastq.clip > $name.fastq.clip.trim &&
	rm $name.fastq $name.fastq.clip

	cat $name.fastq.clip.trim | echo $name.fastq.clip.trim $((`wc -l`/4)) >> $name.temp.log
	gzip $name.fastq.clip.trim
	
	echo "Pre-processing $name... OK" $(date) >> $name.temp.log

}

export -f preprocessing

for i in $files

	do 

	mv $i ./Temp

done
	
cd ./Temp

	ls *.fastq | parallel -j+0 --eta 'preprocessing {}' &&
	cat *.temp.log >> ../log.txt &&
	rm *.temp.log &&
	mv *.fastq.clip.trim.gz ..
		 
cd ..
