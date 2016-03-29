#!/bin/bash
# Edited: 2015/05/19 by Elisabeth D. Chen

# This script reads in raw sequencing data from mulitple samples, processes them separately and sends them into separate pipelines (RNASeq-paired-pipe.sh)
# This script can be executed in a directory containing raw sequencing data. Data from each sample should be in a separate directory.

# Requires: ~/submit.job containing bsub command to submit jobs to farm
#           $TEAM/scripts/RNASeq-paired-pipe.sh with downstream pipeline

echo "Running script..."

for name in $(ls Index* -d) ; do #Read all directories
	#if [ ! -d $filename ] ; then exit ; fi 

	#name=${filename#./}
	echo $name

	cd $name
	ls -lh Alper.$name.?_sequence.txt

	#count=`ls -1 *.fastq.gz 2>/dev/null | wc -l` # Check if unzipped file exists.
	#if [ $count != 0 ]; then gunzip *.fastq.gz ; fi
	#echo "Unzipped $name..."
	#if [ ! -e ${name}_R1.fastq ] ; then cat *L001_R1* *L002_R1* > ${name}_R1.fastq ; fi # Concatanate runs from Lane 1 and Lane 2
	#if [ ! -e ${name}_R2.fastq ] ; then cat *L001_R2* *L002_R2* > ${name}_R2.fastq ; fi
	#echo "Concatanated $name into ${name}_R1.fastq..."
	~/submit.job -n ${name}.process -s ~/my_src/PhD/Alper/RNASeq-paired-pipe_EISA.sh -q long -m 10000 Alper.$name.1_sequence.txt Alper.$name.2_sequence.txt $name 

	cd ..	
done

echo "Finished script."

