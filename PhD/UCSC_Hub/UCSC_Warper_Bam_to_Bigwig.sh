#!/bin/bash


###### Dependencies #####
#						#
#	dos2unix			#
#	fetchchromsizes		#
#	bedToBigBed 		#
# 	beditemoverlapcount	#
#	bedgraphtobigwig	#
#						#
#	hubCheck			#
#	rsync				#
#	Server SSH KEY		#
#						#
#########################


######## INSTALL #######

# conda install --channel https://conda.anaconda.org/trent dos2unix
# conda install --channel https://conda.anaconda.org/bioconda ucsc-fetchchromsizes
# conda install --channel https://conda.anaconda.org/bioconda ucsc-bedtobigbed
# conda install --channel https://conda.anaconda.org/anaconda wget

#conda install --channel https://conda.anaconda.org/bioconda ucsc-beditemoverlapcount
#conda install --channel https://conda.anaconda.org/bioconda ucsc-bedgraphtobigwig



# wget -r http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/hubCheck -O hubCheck && chmod +x hubCheck #&& mv hubCheck $path 
# #wget -r http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/hubCheck -O hubCheck && chmod +x hubCheck #&& mv hubCheck $path 



### Config ###

#All the variables values have to be strings, that is mean that they have to be within " "

Lab="Teichmann"
Project="PV"
Owner="Xiuwei"
email="zhangxiuwei03@gmail.com"


#Only one genome per project is suported
#Use UCSC genome names

Genome="mm10"

#Do not finish path with /
#Use absolute paths

input_dir="/homes/xiuweiz/public_html/bam_files"
output_dir="/homes/xiuweiz/public_html/bam_files/UCSC"
index="/homes/xiuweiz/public_html/bam_files/UCSC/index.txt"

#HTTP server info
#FTP it is slower, do not use it

Server="em-x2.gurdon.cam.ac.uk"
SSH_path="/Library/WebServer/Documents/UCSC"
HTML_path="/UCSC"
User="gp7"
Key="/Users/gp7/emx2_key.pem"



### RUN ###


#Reading index

mkdir -p $Lab 
mkdir -p $Lab/$Owner
mkdir -p $Lab/$Owner/$Project
cd $output_dir/$Lab/$Owner/$Project


####THIS IS THE SH..

#samtools sort -o WT_BAM\@WT_BAM_rep1.sort.bam WT_BAM\@WT_BAM_rep1.bam
#samtools index WT_BAM\@WT_BAM_rep1.sort.bam WT_BAM\@WT_BAM_rep1.sort.bam.bai

while read i;

	do 

		#CSV

		# url=$(echo $i | dos2unix | awk -F "\"*,\"*" '{print $1}')
		# name=$(echo $i | dos2unix| awk -F "\"*,\"*" '{print $2"@"$3"."$4}')
		# Type="$(echo $i | dos2unix| awk -F "\"*,\"*" '{print $4}')" 


		#Tab delimited

		

		url=$(echo $i | dos2unix | awk	'{FS="\t"; print $1}')
		name=$(echo $i | dos2unix| awk	'{FS="\t"; print $2"@"$3"."$4}')
		Type="$(echo $i | dos2unix| awk	'{FS="\t"; print $4}')" 		


		echo $name

		mkdir -p $Type


		if [ -f "./$Type/$name" ]
			then
				echo "${name} found"

			else
				echo "creating $name at $Lab/$Owner/$Project/$Type"
				cp $input_dir/$url $Type/$name
		fi	


	done < $index



fetchChromSizes $Genome > $Genome.chromsizes

echo	hub ${Lab}_Lab_${Owner}	>	hub.txt
echo	shortLabel  ${Lab}/${Owner}/$Project	>>	hub.txt
echo	longLabel ${Lab} Lab - ${Owner} - $Project	>>	hub.txt
echo	genomesFile genomes.txt	>>	hub.txt
echo	email ${email}	>>	hub.txt

cat hub.txt

echo	"genome $Genome"	> genomes.txt
echo	"trackDb trackDb.txt"	>> 	genomes.txt

cat genomes.txt

echo	"#----------------------"	> trackDb.txt
echo 	"# ${Lab} Lab - ${Owner} - ${Project}"	>> trackDb.txt
echo 	"#----------------------"	>> trackDb.txt
echo	>> trackDb.txt
echo	>> trackDb.txt
echo 	"#----------RAW TRACKS----------"	>> trackDb.txt
echo	>> trackDb.txt
echo 	>> trackDb.txt

cat trackDb.txt


#BAMs


Type="bam"
detected=$(ls $Type/)
mkdir -p $Type/Temp
mkdir -p bw


# echo $detected


if [ ! $(ls $Type/ | wc -l) -eq 0 ]

then

	#echo	"$Type detected" $detected	>> trackDb.txt

	echo	>> trackDb.txt
	echo	>> trackDb.txt
	echo 	"#---BAM---"	>> trackDb.txt
	echo	>> trackDb.txt
	echo 	>> trackDb.txt

	echo	track BAM 	>> trackDb.txt
	echo	compositeTrack on 	>> trackDb.txt
	echo	shortLabel BAMs 	>> trackDb.txt
	echo	longLabel BAMs 	>> trackDb.txt
	echo	configurable on 	>> trackDb.txt
	echo	>> trackDb.txt
	echo	>> trackDb.txt
	echo	"### $Type detected###" $(ls $Type)	>> trackDb.txt
	echo	>> trackDb.txt
	echo	>> trackDb.txt




	for i in $(ls $Type/*.$Type)


		do 

			name=$(basename $i .bam)
			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')


			###### -split of bamtobed is the key!!

			echo "@HD	VN:1.4	SO:coordinate" > $Type/$name.head

			awk '{print "@SQ\tSN:"$1"\tLN:"$2}' $Genome.chromsizes >> $Type/$name.head

			#samtools view $Type/$name.bam  > $Type/$name.sam

			samtools view $Type/$name.bam | awk '{OFS="\t"; $3="chr"$3; print}' > $Type/$name.sam


			cat $Type/$name.head $Type/$name.sam | samtools view -Sb - > $Type/$name.bam


			samtools sort -o  $Type/$name.sort.bam $Type/$name.bam
			samtools index $Type/$name.sort.bam $Type/$name.sort.bam.bai


			bamToBed -i $Type/$name.sort.bam -split > $Type/Temp/$name.bed


			sort -k1,1 $Type/Temp/$name.bed | bedItemOverlapCount $Genome -chromSize=$Genome.chromsizes stdin  > $Type/Temp/$name.bedGraph


			awk '$1!="sensor_piRNA_mjIs144"' $Type/Temp/$name.bedGraph > $Type/Temp/$name.bedGraph.filter

			python ~/my_src/PhD/UCSC_Hub/normalize_bw.py $Type/Temp/$name.bedGraph.filter > $Type/Temp/$name.bedGraph.filter.norm

			bedGraphToBigWig $Type/Temp/$name.bedGraph.filter.norm $Genome.chromsizes $Type/Temp/$name.bw



			#bedGraphToBigWig $Type/Temp/$name.bedGraph $Genome.chromsizes $Type/Temp/$name.bw
			
			mv $Type/Temp/$name.bw bw

			#rm $Type/Temp/*
			#rm $Type/*.bam


			# mv $Type/Temp/$name.sort.bam $Type
			# mv $Type/Temp/$name.sorst.bam.bai $Type


			echo	>> trackDb.txt
	        echo	track $name.bam 	>> trackDb.txt
	        echo	parent BAM on 	>> trackDb.txt
	        echo	bigDataUrl $Type/$name.sort.bam 	>> trackDb.txt
	        echo	shortLabel $short_name >> trackDb.txt
	        echo	longLabel $long_name 	>> trackDb.txt
	        echo	type bam 	>> trackDb.txt
	        echo	showNames on	>> trackDb.txt


	 done

	echo	>> trackDb.txt
	echo	>> trackDb.txt

	#rm -rf $Type/Temp/


fi




#BigWigs

Type="bw"
detected=$(ls $Type/)


if [ ! $(ls $Type/ | wc -l) -eq 0 ]

then

	echo	"$Type detected" $detected	

	echo	>> trackDb.txt
	echo	>> trackDb.txt
	echo 	"#---BigWigs---"	>> trackDb.txt
	echo	>> trackDb.txt
	echo 	>> trackDb.txt

	echo	track BW 	>> trackDb.txt
	echo	compositeTrack on 	>> trackDb.txt
	echo	dragAndDrop on 	>> trackDb.txt
	echo	shortLabel BigWigs 	>> trackDb.txt
	echo	longLabel BigWigs 	>> trackDb.txt
	echo	showSubtrackColorOnUi on 	>> trackDb.txt
	echo	configurable on 	>> trackDb.txt
	echo	container multiWig 	>> trackDb.txt
	echo	"type bigWig 0 1000" 	>> trackDb.txt
	echo	>> trackDb.txt
	echo	>> trackDb.txt
	echo	"### $Type detected###" $(ls $Type)	>> trackDb.txt
	echo	>> trackDb.txt
	echo	>> trackDb.txt




	for i in $(ls $Type)


		do 

			name=$(basename $i .bw)
			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')

			#echo $name | grep -o -e'[0-9]\{3\}'
			echo	>> trackDb.txt
	        echo	track $name.bw 	>> trackDb.txt
	        echo	parent BW on 	>> trackDb.txt
	        echo	bigDataUrl $Type/$name.bw 	>> trackDb.txt
	        echo	shortLabel $short_name >> trackDb.txt
	        echo	longLabel $long_name 	>> trackDb.txt
	        echo	type bigWig 	>> trackDb.txt
	        echo	color 0,0,0 	>> trackDb.txt
	        echo	configurable on 	>> trackDb.txt
	        echo	maxHeightPixels 1:100:500 	>> trackDb.txt
	        echo	autoScale on 	>> trackDb.txt


	 done

	echo	>> trackDb.txt
	echo	>> trackDb.txt


fi



##########


# # bed6

# Type="bed6"
# mkdir -p $Type/Temp
# detected=$(ls $Type/)


# if [ ! $(ls $Type/ | wc -l) -eq 0 ]

# 	then

# 		echo	"$Type detected" $files	


# 	for i in $(ls $Type/*.$Type)


# 		do 

# 			name=$(basename $i .$Type)
# 			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
# 			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')

# 			sort -k1,1 -k2,2n $Type/$name.$Type > $Type/Temp/$name.sorted.bed
# 			bedToBigBed $Type/Temp/$name.sorted.bed $Genome.chromsizes $Type/Temp/$name.$Type.bb


# 		done

# fi


# # bed5FloatScore

# Type="bed5FloatScore"
# mkdir -p $Type/Temp
# detected=$(ls $Type/)


# if [ ! $(ls $Type/ | wc -l) -eq 0 ]

# then

# 	echo	"$Type detected" $files	


# 	for i in $(ls $Type/*.$Type)


# 		do 

# 			name=$(basename $i .$Type)
# 			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
# 			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')


# 			awk  '{printf "%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, $6}' $Type/$name.$Type > $Type/Temp/$name.bed
# 			sort -k1,1 -k2,2n $Type/Temp/$name.bed > $Type/Temp/$name.sorted.bed
# 			bedToBigBed $Type/Temp/$name.sorted.bed $Genome.chromsizes $Type/Temp/$name.$Type.bb  

# 		done


# 	# # max_bed_score=$(cat $Type/*.$Type |  awk  'BEGIN{max=0}{if(($5)>max)  max=($5)}END {printf "%.0f\n",  max}')

# fi

# # narrowPeak

# Type="narrowPeak"
# mkdir -p $Type/Temp
# detected=$(ls $Type/)


# if [ ! $(ls $Type/ | wc -l) -eq 0 ]

# then

# 	echo	"$Type detected" $files	


# 	for i in $(ls $Type/*.$Type)


# 		do 

# 			name=$(basename $i .$Type)
# 			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
# 			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')


# 			awk  '{printf "%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $9, $6}' $Type/$name.$Type > $Type/Temp/$name.bed
# 			sort -k1,1 -k2,2n $Type/Temp/$name.bed > $Type/Temp/$name.sorted.bed
# 			bedToBigBed $Type/Temp/$name.sorted.bed $Genome.chromsizes $Type/Temp/$name.$Type.bb  

# 		done


# fi

# # BigBed

# Type="BigBed"
# mkdir -p $Type
# mv	bed6/Temp/*.bb $Type/
# mv  bed5FloatScore/Temp/*.bb $Type/
# mv  narrowPeak/Temp/*.bb $Type/
# detected=$(ls $Type/)




# if [ ! $(ls $Type/ | wc -l) -eq 1 ]

# then

# 	echo This is $(ls $Type/ | wc -l)==0

# 	echo	"$Type detected" $(ls $Type/ | wc -l)

# 	echo	>> trackDb.txt
# 	echo	>> trackDb.txt
# 	echo 	"#---BigBeds---"	>> trackDb.txt
# 	echo	>> trackDb.txt
# 	echo 	>> trackDb.txt


# 	echo	track BB	>> trackDb.txt
# 	echo	compositeTrack on	>> trackDb.txt
# 	echo	longLabel BigBeds	>> trackDb.txt
# 	echo	shortLabel BigBeds	>> trackDb.txt
# 	echo	type bigBed 0 1.0	>> trackDb.txt
# 	echo	viewLimits 0.0:0.2	>> trackDb.txt
# 	echo	visibility pack	>> trackDb.txt
# 	echo	allButtonPair on	>> trackDb.txt
# 	echo	dragAndDrop on	>> trackDb.txt
# 	echo	centerLabelsDense on	>> trackDb.txt
# 	echo	>> trackDb.txt
# 	echo	>> trackDb.txt
# 	echo	"### $Type detected###" $(ls $Type)	>> trackDb.txt
# 	echo	>> trackDb.txt
# 	echo	>> trackDb.txt



# 	for i in $(ls $Type/*.bed6.bb)


# 		do 

# 			name=$(basename $i .bb)
# 			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
# 			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')




# 	        echo	track $name	>> trackDb.txt
# 	        echo	parent BB on	>> trackDb.txt
# 	        echo	bigDataUrl $Type/$name.bb 	>> trackDb.txt
# 	        echo	shortLabel $short_name  	>> trackDb.txt
# 	        echo	longLabel $long_name 	>> trackDb.txt
# 	        echo	type bigBed 6 +	>> trackDb.txt
# 			echo	>> trackDb.txt

# 		done


# 	for i in $(ls $Type/*.bed5FloatScore.bb)


# 		do 

# 			name=$(basename $i .bb)
# 			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
# 			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')

# 	        echo	track $name	>> trackDb.txt
# 	        echo	parent BB on	>> trackDb.txt
# 	        echo	bigDataUrl $Type/$name.bb 	>> trackDb.txt
# 	        echo	shortLabel $short_name  	>> trackDb.txt
# 	        echo	longLabel $long_name 	>> trackDb.txt
# 	        echo	type bigBed 5 .	>> trackDb.txt
# 			echo	>> trackDb.txt

# 		done


# 	for i in $(ls $Type/*.narrowPeak.bb)


# 		do 

# 			name=$(basename $i .bb)
# 			short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
# 			long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')


# 	        echo	track $name	>> trackDb.txt
# 	        echo	parent BB on	>> trackDb.txt
# 	        echo	bigDataUrl $Type/$name.bb 	>> trackDb.txt
# 	        echo	shortLabel $short_name  	>> trackDb.txt
# 	        echo	longLabel $long_name 	>> trackDb.txt
# 	        echo	type bigBed 5 .	>> trackDb.txt
# 			echo	>> trackDb.txt

# 		done

# fi

# cd -




#####################





# rsync_ssh=$(echo ssh -i $Key)

# rsync -Pav -e 'ssh -i /Users/gp7/emx2_key.pem' $Lab $Server:$SSH_path


# hubCheck http://$Server$HTML_path/$Lab/$Owner/$Project/hub.txt

# echo hubCheck http://${Server}${Server_path}/$Lab/$Owner/$Project/hub.txt








# rsync -Pav -e 'ssh -i /Users/gp7/Key.pem.txt' $Lab ubuntu@54.214.245.35:/home/ubuntu/Tracks/


# hubCheck http://54.214.245.35/Tracks/$Lab/$Owner/$Project/hub.txt

# echo hubCheck http://${Server}/${Server_path}/$Lab/$Owner/$Project/hub.txt















####FUTURE####

# # bed5FloatScore

# Type="bed5FloatScore"
# #mkdir -p $Type
# #rsync -aq $input_dir/*.$Type $Type
# detecte=s $Ty/*pe)



# if [ -e "${detected[0]}" ];
# then

# 	echo	"$Type detected" $files	

# echo	>> trackDb.txt
# echo	>> trackDb.txt
# echo 	"#---$Type---"	>> trackDb.txt
# echo	>> trackDb.txt
# echo 	>> trackDb.txt

# echo	track $Type	>> trackDb.txt
# echo	compositeTrack on	>> trackDb.txt
# echo	longLabel BigBeds	>> trackDb.txt
# echo	shortLabel BigBeds	>> trackDb.txt
# echo	type bigBed 0 1.0	>> trackDb.txt
# echo	viewLimits 0.0:0.2	>> trackDb.txt
# echo	visibility pack	>> trackDb.txt
# echo	allButtonPair on	>> trackDb.txt
# echo	dragAndDrop on	>> trackDb.txt
# echo	centerLabelsDense on	>> trackDb.txt
# echo	>> trackDb.txt
# echo	>> trackDb.txt
# echo	"### $Type detected###" $(ls $Type)	>> trackDb.txt
# echo	>> trackDb.txt
# echo	>> trackDb.txt


# for i in $(ls $Type/)


# 	do 

# 		name=$(basename $i .$Type)
		#short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
		#long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')

#         echo	track $name	>> trackDb.txt
#         # echo	parent $Type on	>> trackDb.txt
#         echo	type $Type	>> trackDb.txt        
#         echo	bigDataUrl $Type/$name.$Type	>> trackDb.txt
#         echo	shortLabel $short_name >> trackDb.txt
#         echo	longLabel $long_nameype	>> trackDb.txt
# 		echo	>> trackDb.txt

# 	done


# #narrowPeak

# Type="narrowPeak"
# #mkdir -p $Type
# #rsync -aq $input_dir/*.$Type $Type
# detecte=s $Ty/*pe)



# if [ -e "${detected[0]}" ];
# then

# 	echo	"$Type detected" $files	

# echo	>> trackDb.txt
# echo	>> trackDb.txt
# echo 	"#---$Type---"	>> trackDb.txt
# echo	>> trackDb.txt
# echo 	>> trackDb.txt

# echo	track $Type	>> trackDb.txt
# echo	compositeTrack on	>> trackDb.txt
# echo	longLabel BigBeds	>> trackDb.txt
# echo	shortLabel BigBeds	>> trackDb.txt
# echo	type bigBed 0 1.0	>> trackDb.txt
# echo	viewLimits 0.0:0.2	>> trackDb.txt
# echo	visibility pack	>> trackDb.txt
# echo	allButtonPair on	>> trackDb.txt
# echo	dragAndDrop on	>> trackDb.txt
# echo	centerLabelsDense on	>> trackDb.txt
# echo	>> trackDb.txt
# echo	>> trackDb.txt
# echo	"### $Type detected###" $(ls $Type)	>> trackDb.txt
# echo	>> trackDb.txt
# echo	>> trackDb.txt


# for i in $(ls $Type/)


# 	do 

# 		name=$(basename $i .$Type)
		#short_name=$(echo $name | cut -f1 -d@ | sed 's/_/\ /g')
		#long_name=$(echo $name | cut -f2 -d@ | sed 's/_/\ /g')

#         echo	track $name	>> trackDb.txt
#         # echo	parent $Type on	>> trackDb.txt
#         echo	type $Type	>> trackDb.txt        
#         echo	bigDataUrl $Type/$name.$Type	>> trackDb.txt
#         echo	shortLabel $short_name >> trackDb.txt
#         echo	longLabel $long_nameype	>> trackDb.txt
# 		echo	>> trackDb.txt

# 	done




######ASIA#####


# Lab="Miska"
# Project="piRNA_biogenesis"
# Owner="Asia"
# email="jk571@cam.ac.uk"


# #Only one genome per project is suported
# #Use UCSC genome names

# Genome="ce10"

# #Do not finish path with /
# #Use absolute paths

# input_dir="/Users/gp7/Google_Drive/Results/Asia/piRNA_data_back_up"
# output_dir="/Users/gp7/Google_Drive/USCS_hub"
# index="/Users/gp7/Google_Drive/Results/Asia/piRNA_data_back_up/index.txt"
