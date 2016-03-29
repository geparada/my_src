#!/bin/bash

# Edited on 2015/05/19 by Elisabeth D. Chen
# Edited on 2015/08/03 by Guillermo E. Parada

# This script takes raw sequencing data (one single file per mate) and feeds it through the following process:
# Raw data  >  FastQC  >  STAR aligner (might switch to HISAT)  >  RSeQC to filter rRNA  >  Qualimap  >  HTSeq-count / featureCounts  >  Create summary

# Requires: 1) Raw .fastq files x2 as $1, $2
#           2) Name to create file extensions as $3
#           3) ~/submit.job script to send jobs to the farm
#           4) $TEAM/scripts/extract_data.sh to create final summary

GENOME=$TEAM/Genome/ERCC/ERCC92_ce10/build/ # C. elegans genome index
ECOLI=~/Genomes/ecoli_genome/Ecoli_genome_ENSEMBL 
GTF=$TEAM/Genome/ERCC/ERCC92_ce10/ERCC92_ce10.gtf # C. elegans, annotation file

###CE10

# GENOME=$TEAM/Celegans_genome_ENSEMBL/CE10/build/ # C. elegans genome index
# GTF=$TEAM/Celegans_genome_ENSEMBL/CE10/WS220_gurdon_curated_genes_FC.gtf # C. elegans, annotation file

alignEcoli="No"

while getopts g:e:Eh opt; do
	case $opt in
		g) GENOME=$OPTARG;;
		e) ECOLI=$OPTARG;; # PATH TO ECOLIGENOME
		E) alignEcoli="Yes";; 
		h) echo -e  "This script takes the following commands: \n-g\tspecify C. elegans genome directory. Default is Celegans_genome_ENSEMBL.\n-e\tspecify E. coli genome directory.\n-E\tIf set, alignEcoli=True. Default = False."; exit ;;
		:) echo "Option -$OPTARG requires an argument." >&2 ;;
	esac
done

shift $(( OPTIND -1 )) # Resets command line interface count.

FASTQ1=$1
FASTQ2=$2
NAME=$3

if [ -z $NAME ]; then echo >&2 "Supply a name."; exit ; fi 
if [ -z $FASTQ1 ] ; then echo >&2i $(date) "Enter two valid .fastq files."; exit ; fi
if [ -z $FASTQ2 ] ; then echo >&2i $(date) "Enter two valid .fastq files."; exit ; fi



###### Pre-procesing #########

echo "Pre-processing"

raw_count=$(cat $FASTQ1 | echo $((`wc -l`/4)))

if [ ! -d "Pre-processing" ]; then mkdir Pre-processing; fi

cutadapt -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -m 35 -o Pre-processing/$FASTQ1 -p Pre-processing/$FASTQ2 $FASTQ1 $FASTQ2 > $NAME.cutadapt.summary

cd Pre-processing

	sickle pe -f $FASTQ1 -r $FASTQ2 -t sanger -q 20 -l 35 -o ../$FASTQ1 -p ../$FASTQ2 -s unpaired.fastq
cd ..

pre_processed_count=$(cat $FASTQ1 | echo $((`wc -l`/4)))
survivior_reads=$(bc -l <<< "$pre_processed_count*100/$raw_count")

echo "Raw reads = $raw_count" > $NAME.pre-processing.summary
echo "Reads after pre-processing = $pre_processed_count (% $survivior_reads)" >> $NAME.pre-processing.summary

#                     FastQC
echo $(date) "Running FastQC"
# Run FastQC and store results in directory FastQC.results
if [ ! -d "FastQC.results" ]; then mkdir FastQC.results; fi

~/submit.job -n $NAME.fastqc -s $TEAM/scripts/.fastqc.sh -m 1000 $FASTQ1 $FASTQ2


#                     STAR

# Add STAR to the path
#export PATH=$PATH:/nfs/users/nfs_e/ec12/Programs/STAR_2.4.0j/source/


if [ $alignEcoli == "Yes" ]; then 
	echo $(date) "Aligning to E coli genome";
	# Check if STAR.ecoli folder already exists. If not, make the folder. Change directory into the folder
	if [ ! -d "STAR.ecoli.results" ]; then mkdir STAR.ecoli.results ; fi
	cd STAR.ecoli.results ;
	# Code to run STAR alignment. Use ThreadN 28. Also create new fastq file with unmapped reads
	STAR --genomeDir $ECOLI --readFilesIn ../$FASTQ1 ../$FASTQ2 --runThreadN 28 --outFileNamePrefix $NAME. --outReadsUnmapped Fastx ;

	cd .. ; fi


echo $(date) "Aligning to C elegans genome"
# Check if STAR.worm folder already exists. If not, make folder. Change directory into the folder
if [ ! -d "STAR.worm.results" ]; then mkdir STAR.worm.results ; fi 

cd STAR.worm.results

	if [ $alignEcoli == "Yes" ]; then 
		FASTQ1=STAR.ecoli.results/$NAME.Unmapped.out.mate1
		FASTQ2=STAR.ecoli.results/$NAME.Unmapped.out.mate2
	fi


	# Code to run STAR alignment. 
	STAR --genomeDir $GENOME --readFilesIn ../$FASTQ1 ../$FASTQ2 --runThreadN 28 --outFileNamePrefix $NAME. --outReadsUnmapped Fastx  --outSAMattributes All

	awk '($18 ~ 0) { print > "non-canonical"; next } { print > "canonical" }' $NAME.Aligned.out.sam
	mv non-canonical $NAME.Aligned.out.sam.NC_SJ
	mv canonical $NAME.Aligned.out.sam.C_SJ

	SAM=$NAME.Aligned.out.sam.C_SJ

	if [ ! -e $SAM ]; then echo >&2 "Alignment failed $SAM is not a file."; exit ; fi

cd ..

#                      QUALIMAP & RSeQC
echo $(date) "Sorting SAM..."
samtools view -bS STAR.worm.results/$SAM | samtools sort -n - $SAM.sorted

BAM=$SAM.sorted.bam

if [ ! -e $BAM ]; then echo >&2 "Sorting did not work. $BAM is not a file." ; exit ; fi

echo $(date) "Submitted qualimap job"

~/submit.job -n qualimap$NAME -s $TEAM/scripts/.qualimap.sh -m 2000 $BAM $GTF
echo $(date) "Filtering rRNA..."
$TEAM/scripts/.RSeQC.sh $BAM $NAME


#         featureCounts

 ~/Programs/subread-1.4.6-p2-source/bin/featureCounts -p -T 5 -a $TEAM/Celegans_genome_ENSEMBL/Celegans_alignment_ENSEMBL.gtf -o ${NAME}.featureCounts STAR.worm.results/$SAM > ${NAME}.featureCounts.interface

awk 'NR>2{ print $1, $7 }' OFS="\t" $NAME.featureCounts > $NAME.count.out


#         HTSEQ-COUNT

BAM=$NAME.ex.bam

if [ ! -e RSeQC.results/$BAM ]; then echo >&2 "Filtering rRNA did not work. $BAM is not a file." ; exit ; fi

echo $(date) "Counting reads with HTSeq-count..."
samtools view RSeQC.results/$BAM | ~/Programs/HTSeq-0.6.1p1/scripts/htseq-count --order=name - $GTF > ./$NAME.htseq-count.out

$TEAM/scripts/extract_data.sh $NAME

echo -e $(date) "\n\nPipeline finished."


echo $NAME | mutt -s "DONE" -a PE40.exons.featureCounts -- geparada88@gmail.com
