Genome="/home/geparada/OLego/hg19"
TAGs="/media/HD3/Resultados/Micro_exons/Tags/SJ_tags_exon_25_micro_exons.fa"
#TAGs="/media/HD3/Resultados/Micro_exons/Tags/SJ_tags_exon_25_micro_exons.canonical.fa"
TAGs_round2="/media/HD3/Resultados/Micro_exons/Tags/Round2/ME_canonical_SJ_tags.fa"
READS_FILES="/media/HD4/db/Illumina_bodymap/*trim.gz"
#READS_FILES="/media/HD4/db/Illumina_bodymap/Tissues/*.trim.gz"
#READS_FILES="/media/HD4/db/si_UPF1/*.gz"


#############################################################################

                     #De novo Olego/processing

#############################################################################



# cd De_novo
# #cd Test

# for i in $READS_FILES

# 	do

# 	gzip -d $i                                          #Necesario para los que estan comprimidos



# 	#name=$(echo $i | xargs -i basename {} )
# 	name=$(echo ${i%.*} | xargs -i basename {} )        #Necesario para los que estan comprimidos

# 	echo $name
# 	#olego ~/OLego/hg19 -t 8 $i > $name.sam
# 	olego ~/OLego/hg19 -t 8 ${i%.*} > $name.sam	        #Necesario para los que estan comprimidos



# 	samtools view -SH $name.sam | grep "|" -v > $name.sam.uniq.SJ   #header

# 	samtools view -S $name.sam | awk 'NF==20 || NF==21'| awk '$15=="X0:i:1"' > $name.sam.uniq
# 	awk '($6 ~ /N/)' $name.sam.uniq >> $name.sam.uniq.SJ
# 	samtools view -Sb $name.sam.uniq.SJ -o $name.bam.uniq.SJ


# 	number_of_mapped_reads=$(samtools view -S $name.sam | awk 'NF==20 || NF==21' | wc -l)
# 	number_of_unmapped_reads=$(samtools view -S $name.sam | awk '$2==4' | wc -l )
# 	number_of_uniq_mapped_reads=$(cat $name.sam.uniq | wc -l)
# 	number_of_uniq_SJ_mapped_reads=$(samtools view -S $name.sam.uniq.SJ | wc -l)

# 	echo "Number of mapped reads:" $number_of_mapped_reads
# 	echo "Number of unmapped reads:" $number_of_unmapped_reads
# 	echo "Number of uniq mapped reads:" $number_of_uniq_mapped_reads
# 	echo "Number of uniq SJ mapped reads:" $number_of_uniq_SJ_mapped_reads

	
# 	samtools view -Sb $name.sam| samtools sort - $name
# 	bamToBed -i $name.bam -split > $name.bed

	
# 	awk '{if($6=="+") print}' $name.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n > $name.plus.bedGraph &&
# 	awk '{if($6=="-") print}' $name.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > $name.minus.bedGraph &&


# 	bedGraphToBigWig $name.plus.bedGraph chromInfo.txt $name.plus.bw &&
# 	bedGraphToBigWig $name.minus.bedGraph chromInfo.txt $name.minus.bw &&


# 	rm $name.sam.uniq $name.bam *bed* $name.sam.uniq.SJ $name.sam

# 	gzip ${i%.*}                                         #Necesario para los que estan comprimidos

# 	done


# cd ..


#############################################################################

                     #Round1 - Tag-based aligment and pre-processing

#############################################################################






#cd Tags/Round1/IBM/Mixture
cd Tags/Round1/IBM/Tissue
#cd Tags/Round1/shUPF1/

ls $READS_FILES

 for i in $(ls $READS_FILES)

 	do

	####MAPPING#####

	name=$(echo ${i%%.*} | xargs -i basename {} )        #Necesario para los que estan comprimidos

	fastq-dump $name.sra

	# zcat $i > $name.fastq &&                                         #Necesario para los que estan comprimidos
	#ln -s $i $name.fastq

	echo $name

	bwa mem -t 8 -O 2,6 -L 25 $TAGs $name.fastq  > $name.sam
	python ~/my_src/ME/Pipeline/Round_1/alingment_pre_processing.py $name.sam F > $name.sam.pre_processed &&

	echo "OK"

	rm $name.sam $name.fastq


# 	# ###POST-Processing###

# 	python ~/my_src/ME/Pipeline/Round_1/row_ME.py ~/db/genome/hg19.fa $name.sam.pre_processed > $name.sam.row_ME
# 	python ~/MapSplice-v2.1.8/mapsplice.py -c ~/db/genome/hg19 -x ~/MapSplice-v2.1.8/bowtie_index_hg19 -1 $name.sam.row_ME.fastq -p 8
# 	mv mapsplice_out/alignments.sam $name.sam.row_ME.hg19.sam
# 	python ~/my_src/Tools/FASTQ_to_FASTA.py $name.sam.row_ME.fastq > $name.sam.row_ME.fasta

# 	dust $name.sam.row_ME.fasta | grep NNNNN -B 1 | grep ^\>  > $name.sam.row_ME.fastq.dust
# 	blat ~/RepeatMasker/Libraries/Homo_sapiens_all_Repbase $name.sam.row_ME.fasta $name.sam.row_ME.fastq.Repbase -fastMap -noHead

# 	python ~/my_src/ME/Pipeline/Round_1/ME_filter1.py ~/db/genome/hg19.fa $name.sam.row_ME $name.sam.row_ME.hg19.sam $name.sam.row_ME.fastq.dust $name.sam.row_ME.fastq.Repbase ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/hg19.100way.phyloP100way.bw ~/db/hg19.46way.phyloP46way.primates.bw > $name.sam.row_ME.filter1

# 	done

# #cd ../../../..
# cd ../../..


#############################################################################

                     #Round1 - Micro exon filters 

#############################################################################


# cd Tags/Round1

# #### cat IBM/Mixture/*.filter1 IBM/Tissue/*filter1 shUPF1/*filter1 > TOTAL.sam.row_ME.filter1
# cat IBM/Mixture/*.filter1 IBM/Tissue/*filter1 > TOTAL.sam.row_ME.filter1
# python ~/my_src/ME/Pipeline/Round_1/ME_centric_table.py TOTAL.sam.row_ME.filter1 > TOTAL.sam.row_ME.filter1.ME_centric
# python ~/my_src/ME/Pipeline/Round_1/micro_exons_gencode_row.py ~/db/genome/hg19.fa ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12  ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/hg19.100way.phyloP100way.bw ~/db/hg19.46way.phyloP46way.primates.bw TOTAL.sam.row_ME.filter1.ME_centric > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode
# python ~/my_src/ME/Pipeline/Round_1/ME_filter2.py ~/db/genome/hg19.fa TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode ~/db/Variation/snp138Common.fix ../../../Non_canonical_introns/TOTAL/non_canonical | sort -k 1 -n -r > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2
# python ~/my_src/ME/Pipeline/Round_1/ME_filter3.py TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2 ../../simulation/ROC/ME_sim.scores ~/db/transcriptome/hg19/Gene_models/gencode/v17/gencode.v17.annotation.bed12 > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3

# cd ../..

# #############################################################################

#                      #Round2

# #############################################################################



# cd Tags/Round2

# python ~/my_src/ME/Pipeline/Round_2/Micro_exons_tags.py $TAGs ../Round1/TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3 > TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3.ME_tags.fa
# cat TOTAL.sam.row_ME.filter1.ME_centric.MEs_gencode.filter2.filter3.ME_tags.fa $TAGs> ME_canonical_SJ_tags.fa
# bwa index ME_canonical_SJ_tags.fa

# cd ../..

#####MAPPING#####



#cd Tags/Round2/shUPF1/
# cd Tags/Round2/IBM/Mixture

# for i in $READS_FILES

#  	do

# 	name=$(echo ${i%%.*} | xargs -i basename {} )        #Necesario para los que estan comprimidos
# 	echo $name



# 	####### Maping ######

# 	zcat $i > $name.fastq &&
# 	bwa mem -t 8 -O 2,6 -L 25 $TAGs_round2 $name.fastq  > $name.sam &&
# 	python ~/my_src/ME/Pipeline/Round_2/alingment_pre_processing_round2.py $name.sam F > $name.sam.pre_processed
# 	rm $name.sam $name.fastq

# # 	####POST-Processing###

# 	python ~/my_src/ME/Pipeline/Round_2/round2_ME_reads_fastq.py $name.sam.pre_processed > $name.sam.pre_processed.fastq
# 	bowtie ~/db/genome/index/bowtie1/bowtie_index_hg19 -p 7 -q $name.sam.pre_processed.fastq -S | awk '$2==0 || $2==16'> $name.sam.pre_processed.hg19.sam
# 	python ~/my_src/Tools/FASTQ_to_FASTA.py $name.sam.pre_processed.fastq > $name.sam.pre_processed.fasta
# 	dust $name.sam.pre_processed.fasta | grep NNNNN -B 1 | grep ^\>  > $name.sam.pre_processed.fastq.dust
# 	blat ~/RepeatMasker/Libraries/Homo_sapiens_all_Repbase $name.sam.pre_processed.fasta $name.sam.pre_processed.fastq.Repbase -fastMap -noHead
# 	python ~/my_src/ME/Pipeline/Round_2/Filter1_round2.py $name.sam.pre_processed $name.sam.pre_processed.fastq.dust $name.sam.pre_processed.fastq.Repbase $name.sam.pre_processed.hg19.sam > $name.sam.pre_processed.filter1
	

# # #rm $name.fastq

# done

# #cd ../../..
# cd ../../../..


##############################

		#Final  ---- 	ARRERGLAR TAGS Y ME_SJ_COVERAGE! QUE LOS TAGS DIGAN DE QUE ME SON!

##############################

# cd Tags/Round2/IBM/Tissue

# for i in $(ls *1x75*pre_processed)
# 	do 
# 	name=${i%%1x75*};echo $name
# 	cat $name\1x75.sam.pre_processed.filter1 $name\2x50.sam.pre_processed.filter1 > TOTAL/$name.TOTAL.sam.pre_processed.filter1
# 	done

# cd TOTAL

# for i in $(ls *.TOTAL.sam.pre_processed.filter1)
# 	do 
# 	python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 $i > $i.ME_SJ_coverage
# 	rm $i
# 	done

# cd ../

# cd ../../../..

# #######################----------------------

# cd Tags/Round2/shUPF1/

# cat HepG2_control*.filter1 > TOTAL/HepG2_control.TOTAL.sam.pre_processed.filter1
# cat HepG2_shUPF2_*.filter1 > TOTAL/HepG2_shUPF2.TOTAL.sam.pre_processed.filter1

# cat sicontrol_*.filter1 > TOTAL/sicontrol.TOTAL.sam.pre_processed.filter1
# cat siUPF*.filter1 > TOTAL/siUPF1.TOTAL.sam.pre_processed.filter1

# cat HepG2_control_TOTAL.filter1 sicontrol_TOTAL.filter1 > TOTAL/control.TOTAL.sam.pre_processed.filter1
# cat HepG2_shUPF2_TOTAL.filter1 siUPF1_TOTAL.filter1 > TOTAL/sh.TOTAL.sam.pre_processed.filter1

# cd TOTAL

# for i in $(ls *.TOTAL.sam.pre_processed.filter1)
# 	do 
# 	python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round2/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3.ME_tags.fa /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 $i > $i.ME_SJ_coverage
# 	rm $i
# 	done

# cd ..

# cd ../../..

######################----------------------

cd Tags/Round2/IBM/Mixture

#cat ERR0308*.filter1 > Mixture.TOTAL.filter1
python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 Mixture.TOTAL.filter1 > Mixture.TOTAL.filter1.ME_SJ_coverage

cd ../../..






# #python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round2/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3.ME_tags.fa /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 ~/db/transcriptome/hg19/Gene_models/gencode/v19/gencode.v19.annotation.bed12 HepG2_control_TOTAL.filter1 > HepG2_control_TOTAL.filter1.ME_SJ_coverage

# for i in $(ls *TOTAL*);do echo $i; python ~/my_src/ME/Pipeline/Round_2/ME_SJ_coverage.py /media/HD3/Resultados/Micro_exons/Tags/Round2/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3.ME_tags.fa /media/HD3/Resultados/Micro_exons/Tags/Round1/TOTAL.sam.row_ME.filter1.ME_centric.filter2.filter3 $i> $i.Coverage_Table ;done



# #cat HepG2_control*.filter1 > HepG2_control_TOTAL.filter1

# cd ../../..

# cd Tags/Round2/Final

# cat ../IBM/Mixture/*.filter1 ../IBM/Tissue/*.filter1 ../shUPF1/*.filter1 > TOTAL.Round2.filter1

