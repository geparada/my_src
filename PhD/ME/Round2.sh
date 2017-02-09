#############################################################################

                     #Round2

#############################################################################


TAGs="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/ME/Round1/mm10.ME_TAGs.fa"
TAGs_round2="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/ME/Round2/mm10.ME_TAGs.round1-2.fa"

Genome_index="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/mm10_STAR_build"


Genome_index_bowtie="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/mm10_bowtie_index/mm10" 

Genome="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/mm10.fa"

Genome_name="mm10"


GT_AG_U2_5="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix"
GT_AG_U2_3="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix"

vertebrates_phylop="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60way.bw"
close_phylop="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60wayPlacental.bw"

Gene_anontation_bed12="/lustre/scratch108/compgen/team218/gp7/Genome/hg19/Tracks/Gene_annotation/gencode.v19.chr_patch_hapl_scaff.annotation.bed12"


fastq="/lustre/scratch108/compgen/team218/gp7/Micro-exons/Mouse_ENCODE/fastq"


Round1="/lustre/scratch108/compgen/team218/gp7/Micro-exons/Mouse_ENCODE/Round1"

# mkdir -p Round2


cd /lustre/scratch108/compgen/team218/gp7/Micro-exons/Mouse_ENCODE/Round2

# python ~/my_src/ME/Pipeline/Round_2/Micro_exons_tags.py $TAGs ../Round1/TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3 > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3.ME_tags.fa
# cat TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3.ME_tags.fa $TAGs> ME_canonical_SJ_tags.fa
# bwa index ME_canonical_SJ_tags.fa

# cd ../..

# ####MAPPING#####



# cd Tags/Round2/shUPF1/
# cd Tags/Round2/IBM/Mixture

# for i in $READS_FILES

#  	do

# 	name=$(echo ${i%%.*} | xargs -i basename {} )        #Necesario para los que estan comprimidos
# 	echo $name



# 	####### Maping ######

# 	zcat $fastq/NAME.fastq.gz > NAME.fastq
# 	#bwa mem -t 5 -O 2,6 -L 25 $TAGs_round2 NAME.fastq  > NAME.sam

# 	bowtie $TAGs_round2 -p 5 -q NAME.fastq -S -v 2 > NAME.sam


# 	rm NAME.fastq

# 	# python ~/my_src/ME/Pipeline/Round_2/alingment_pre_processing_round2.py NAME.sam STRANDED > NAME.sam.pre_processed

# 	python ~/my_src/ME/Pipeline/Round_2/alingment_pre_processing_round2_bowtie.py NAME.sam STRANDED > NAME.sam.pre_processed



# 	rm NAME.sam

# # 	####POST-Processing###

#  	python ~/my_src/ME/Pipeline/Round_2/round2_ME_reads_fastq.py NAME.sam.pre_processed > NAME.sam.pre_processed.fastq
#  	bowtie $Genome_index_bowtie -p 5 -q NAME.sam.pre_processed.fastq -S -v 2| awk '$2==0 || $2==16'> NAME.sam.pre_processed.hg19.sam

#  	rm NAME.sam.pre_processed.fastq

# # 	python ~/my_src/Tools/FASTQ_to_FASTA.py NAME.sam.pre_processed.fastq > NAME.sam.pre_processed.fasta
# # 	dust NAME.sam.pre_processed.fasta | grep NNNNN -B 1 | grep ^\>  > NAME.sam.pre_processed.fastq.dust
# # 	blat ~/RepeatMasker/Libraries/Homo_sapiens_all_Repbase NAME.sam.pre_processed.fasta NAME.sam.pre_processed.fastq.Repbase -fastMap -noHead
#  	python ~/my_src/ME/Pipeline/Round_2/Filter1_round2.py NAME.sam.pre_processed NAME.sam.pre_processed.hg19.sam > NAME.sam.pre_processed.filter1

#  	rm NAME.sam.pre_processed NAME.sam.pre_processed.hg19.sam
	

	python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py $TAGs_round2 $Round1/TOTAL.sam.row_ME.filter1.ME_centric $Gene_anontation_bed12 NAME.sam.pre_processed.filter1 NAME > NAME.filter1.ME_SJ_coverage





# done

# #cd ../../..
# cd ../../../..
