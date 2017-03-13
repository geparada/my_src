
####HUMAN####


# TAGs="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/hg19.ME_TAGs.fa"

# Genome_index="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/hg19_STAR_build"

# Genome="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/hg19.fa"

# Genome_name="hg19"

# RepBase="/lustre/scratch117/cellgen/team218/gp7/Genome/RepBase21.08.fasta/humrep_sub.ref"

# hg19_GT_AG_U2_5="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/Tracks/SpliceRack/hg19_GT_AG_U2_5.good.matrix"
# hg19_GT_AG_U2_3="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/Tracks/SpliceRack/hg19_GT_AG_U2_3.good.matrix"

# vertebrates_phylop="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/Tracks/Phylop/hg19.100way.phyloP100way.bw"
# primates_phylop="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/Tracks/Phylop/hg19.46way.phyloP46way.primates.bw"

# Gene_anontation_bed12="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/Tracks/Gene_annotation/gencode.v19.chr_patch_hapl_scaff.annotation.bed12"




####MOUSE####



TAGs="/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/ME/Round1/mm10.ME_TAGs.fa"

Genome_index="/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/mm10_STAR_build"

Genome="/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/mm10.fa"

Genome_name="mm10"


GT_AG_U2_5="/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix"
GT_AG_U2_3="/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix"

vertebrates_phylop="/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60way.bw"
close_phylop="/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60wayPlacental.bw"

Gene_anontation_bed12="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/Tracks/Gene_annotation/gencode.v19.chr_patch_hapl_scaff.annotation.bed12"
Gene_anontation_GTF="/lustre/scratch117/cellgen/team218/gp7/Genome/hg19/Tracks/Gene_annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"


#############################################################################

                     #MAPPING#####

#############################################################################


####MAPPING#####


echo NAME

# fastq-dump NAME.sra

zcat NAME.fastq.gz >  NAME.fastq

# #bwa mem -t 5 -O 2,6 -L 25 $TAGs NAME.fastq  > NAME.sam
# #python ~/my_src/ME/Pipeline/Round_1/alingment_pre_processing.py NAME.sam STRANDED > NAME.sam.pre_processed 
# # echo "OK"

# #rm NAME.sam


mkdir -p Gene_count



#STAR --genomeDir $Genome_index --readFilesIn NAME.fastq --runThreadN 5 --outFileNamePrefix  Gene_count/NAME.sam


cd Gene_count

python ~/my_src/PhD/Teichmann/Get_introns_from_sam.py NAME.samAligned.out.sam 8 > NAME.introns

# featureCounts -a $Gene_anontation_GTF -o NAME.featurecounts.gene NAME.samAligned.out.sam -T 5

# featureCounts -a $Gene_anontation_GTF -o NAME.featurecounts.exon NAME.samAligned.out.sam -T 5 -f -J -G $Genome




# rm NAME.sam

cd ..

rm NAME.fastq







#############################################################################

                     #POST-Processing###

#############################################################################


# mkdir -p Round1

# cd Round1


# python ~/my_src/ME/Pipeline/Round_1/row_ME.py $Genome ../NAME.sam.pre_processed > NAME.sam.row_ME


# #STAR --genomeDir $Genome_index --readFilesIn NAME.sam.row_ME.fastq --runThreadN 1 --outFileNamePrefix NAME.sam.row_ME.$Genome_name. 



# # #python ~/MapSplice-v2.1.8/mapsplice.py -c ~/db/genome/hg19 -x ~/MapSplice-v2.1.8/bowtie_index_hg19 -1 NAME.sam.row_ME.fastq -p 8
# # # mv mapsplice_out/alignments.sam NAME.sam.row_ME.hg19.sam

# # python ~/my_src/Tools/FASTQ_to_FASTA.py NAME.sam.row_ME.fastq > NAME.sam.row_ME.fasta

# mv NAME.sam.row_ME.$Genome_name.Aligned.out.sam NAME.sam.row_ME.$Genome_name.sam

# # dustmasker -in NAME.sam.row_ME.fasta -out NAME.sam.row_ME.dustmasker.fasta
# # grep -e ^'[0-9]' -B 1 NAME.sam.row_ME.dustmasker.fasta | grep ^\> > NAME.sam.row_ME.dustmasker.fasta.read_IDs

# #dust NAME.sam.row_ME.fasta | grep NNNNN -B 1 | grep ^\>  > NAME.sam.row_ME.fastq.dust
# #blat $RepBase NAME.sam.row_ME.fastq.Repbase -fastMap -noHead

# ##python ~/my_src/ME/Pipeline/Round_1/ME_filter1.py ~/db/genome/hg19.fa NAME.sam.row_ME NAME.sam.row_ME.$Genome_name.sam NAME.sam.row_ME.dustmasker.fasta.read_IDs NAME.sam.row_ME.fastq.Repbase ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/hg19.100way.phyloP100way.bw ~/db/hg19.46way.phyloP46way.primates.bw > NAME.sam.row_ME.filter1

# python ~/my_src/ME/Pipeline/Round_1/ME_filter1.py $Genome NAME.sam.row_ME NAME.sam.row_ME.$Genome_name.sam  $GT_AG_U2_5 $GT_AG_U2_3 $vertebrates_phylop $close_phylop> NAME.sam.row_ME.filter1


#############################################################################

                     #Round1 - Micro exon filters 

#############################################################################



# mkdir TOTAL

# cat *.sam.row_ME.filter1 > TOTAL/TOTAL.sam.row_ME.filter1 

# cd TOTAL

# 	python ~/my_src/ME/Pipeline/Round_1/ME_centric_table.py TOTAL.sam.row_ME.filter1 > TOTAL.sam.row_ME.filter1.ME_centric
#	python ~/my_src/ME/Pipeline/Round_2/Micro_exons_tags.py $TAGs TOTAL.sam.row_ME.filter1.ME_centric  > TOTAL.sam.row_ME.filter1.ME_centric.ME_tags.fa

#	mkdir -p ../../Round2

#	cat TOTAL.sam.row_ME.filter1.ME_centric.ME_tags.fa $TAGs> ../../Round2/ME_canonical_SJ_tags.fa
#	cd ../../Round2
#	bowtie-build ME_canonical_SJ_tags.fa ME_canonical_SJ_tags.fa




##### NOT IN USE ######

# 	python ~/my_src/ME/Pipeline/Round_1/micro_exons_gencode_row.py $Genome $Gene_anontation_bed12  $hg19_GT_AG_U2_5 $hg19_GT_AG_U2_3 $vertebrates_phylop $primates_phylop TOTAL.sam.row_ME.filter1.ME_centric > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode
# 	python ~/my_src/ME/Pipeline/Round_1/ME_filter2.py ~/db/genome/hg19.fa TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode ~/db/Variation/snp138Common.fix ../../../Non_canonical_introns/TOTAL/non_canonical | sort -k 1 -n -r > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2
# 	python ~/my_src/ME/Pipeline/Round_1/ME_filter3.py TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2 ../../simulation/ROC/ME_sim.scores ~/db/transcriptome/hg19/Gene_models/gencode/v17/gencode.v17.annotation.bed12 > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3

#####################


# cd ../..
