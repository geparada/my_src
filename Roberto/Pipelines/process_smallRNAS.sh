#for i in $(ls *.fastq)

files=$(echo wgEncodeCshlShortRnaSeqA549CytosolShorttotalCiptapRawDataRep3.fastq wgEncodeCshlShortRnaSeqA549CytosolShorttotalTapRawDataRep3.fastq wgEncodeCshlShortRnaSeqA549NucleusShorttotalCiptapRawDataRep3.fastq wgEncodeCshlShortRnaSeqA549NucleusShorttotalTapRawDataRep3.fastq)


for i in $files

do
name=${i%.fastq}

#echo "Pre-processing..." $(date)
#cat $name.fastq | echo $name.fastq $((`wc -l`/4))
#fastx_clipper -a AAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -l 22 -i $name.fastq -o $name.fastq.clip || fastx_clipper -a AAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -l 22 -i $name.fastq -o $name.fastq.clip -Q 33
#python ~/my_src/Roberto/post_clip_trim_smallRNAs.py $name.fastq.clip > $name.fastq.clip.trim
#cat $name.fastq.clip | echo $name.fastq.clip $((`wc -l`/4))
#cat $name.fastq.clip.trim | echo $name.fastq.clip.trim $((`wc -l`/4))

echo "Processing $name..." $(date)

gzip -d $name.fastq.clip.trim.gz
python ~/my_src/Roberto/CG_adapter.py $name.fastq.clip.trim > $name.fastq.clip.trim.CG &&
cat $name.fastq.clip.trim.CG | echo $name.fastq.clip.trim.CG $((`wc -l`/4)) &&
gzip $name.fastq.clip.trim

echo "Mapping ..." $(date)
bowtie tags_and_hg19 --trim5 6 -a --best --strata -v2 -M10 -S -p8 -q $name.fastq.clip.trim.CG $name.fastq.clip.trim.CG.sam 
echo "Mapping without trimming the adapter ..." $(date)
bowtie tags_and_hg19 -a --best --strata -v2 -M10 -S -p8 -q $name.fastq.clip.trim.CG | awk '$2==0||$2==16' > $name.fastq.clip.trim.CG.sam.adapter && rm $name.fastq.clip $name.fastq.clip.trim.CG


echo "Post-processing..." $(date)
samtools view -SH $name.fastq.clip.trim.CG.sam | grep "|" -v > $name.fastq.clip.trim.CG.sam.header
python ~/my_src/Roberto/process_SAM_smallRNAs.py ~/db/genome/hg19.fa $name.fastq.clip.trim.CG.sam $name.fastq.clip.trim.CG.sam.adapter Galaxy6-\[Filter_sequences_by_length_on_data_5\].fasta ~/db/transcriptome/hg19/Gene_models/RefSeq/refMrna.fa ~/db/transcriptome/hg19/Gene_models/RefSeq/RefSeq_gen_names ~/db/transcriptome/hg19/Gene_models/RefSeq/RefSeqLinks ~/db/transcriptome/hg19/Gene_models/RefSeq/RefSeq.bed12 &&
##Outputs =  temporal_best.sam, antisense_SJ.sam, antisense_SJ_uniq.sam, sense_SJ.sam, sense_SJ_uniq.sam###
rm $name.fastq.clip.trim.CG.sam


	for out in $(echo antisense_SJ.sam antisense_SJ_uniq.sam sense_SJ.sam sense_SJ_uniq.sam)

	do


		with_genomic_adapter=$(awk '$(NF)=="XM:Z:T"' $out | wc -l)
		with_poliA3=$(awk '$(NF-1)=="XP:Z:T"' $out | wc -l)				

		wc -l $out
		echo "Whith a genomic '5 CG barcode:" $with_genomic_adapter
		echo "With a 3' poliA at mRNA:" $with_poliA3
		
		awk '$(NF)=="XM:Z:F" && $(NF-1)=="XP:Z:F"' $out > $out.filter
		
		cat $out.filter | echo "Final:" $((`wc -l`))
		
		awk '{print length($10)}' $out.filter  | textHistogram -maxBinCount=100 stdin
		
		cp $out ./SJ/$name.TOTAL.$out
	
	
	done


###Generando BAMs de alineamientos a SJ separados por sense/antisense total/uniq###

cat $name.fastq.clip.trim.CG.sam.header antisense_SJ.sam.filter | samtools view -Sb - -o temporal.bam && samtools sort temporal.bam $name.antisense_SJ && samtools index $name.antisense_SJ.bam $name.antisense_SJ.bam.bai
cat $name.fastq.clip.trim.CG.sam.header antisense_SJ_uniq.sam.filter | samtools view -Sb - -o temporal.bam && samtools sort temporal.bam $name.antisense_SJ_uniq && samtools index $name.antisense_SJ_uniq.bam $name.antisense_SJ_uniq.bam.bai
cat $name.fastq.clip.trim.CG.sam.header sense_SJ.sam.filter | samtools view -Sb -  -o temporal.bam && samtools sort temporal.bam $name.sense_SJ && samtools index $name.sense_SJ.bam $name.sense_SJ.bam.bai
cat $name.fastq.clip.trim.CG.sam.header sense_SJ_uniq.sam.filter | samtools view -Sb - -o temporal.bam && samtools sort temporal.bam $name.sense_SJ_uniq && samtools index $name.sense_SJ_uniq.bam $name.sense_SJ_uniq.bam.bai

###Generando BAMs de todos los alineamientos separados por total/uniq###
grep "|chr" -v temporal_best.sam | awk '$(NF-2)=="XB:Z:T" && $(NF)=="XM:Z:F" && $(NF-1)=="XP:Z:F"' | cat $name.fastq.clip.trim.CG.sam.header *SJ.sam.filter - | samtools view -Sb - -o $name.TOTAL.all.bam
grep "|chr" -v temporal_best.sam | awk '$(NF-2)=="XB:Z:T" && $(NF)=="XM:Z:F" && $(NF-1)=="XP:Z:F"' | awk '$5!=0' | rev | uniq -u -f 13 | rev | cat $name.fastq.clip.trim.CG.sam.header *SJ_uniq.sam.filter - | samtools view -Sb - -o $name.TOTAL.uniq.bam  

###Separando en carpetas###
mv $name.*sense*.bam* ./SJ
mv $name.TOTAL*.bam ./TOTAL

rm $name.fastq.clip.trim.CG.sam.adapter $name.fastq.clip.trim.CG.sam.header

cd SJ

##Generando BigWigs de alineamientos a SJ###
samtools merge $name.SJ.all.bam $name.antisense_SJ.bam $name.sense_SJ.bam
samtools merge $name.SJ.uniq.bam $name.antisense_SJ_uniq.bam $name.sense_SJ_uniq.bam

	for i in $(ls $name.SJ*.bam)

	do

	name2=${i%.bam}
	echo $name2

	bamToBed -i $name2.bam -split > $name2.bed

	awk '{if($6=="+") print}' $name2.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n > $name2.plus.bedGraph
	awk '{if($6=="-") print}' $name2.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > $name2.minus.bedGraph

	bedGraphToBigWig $name2.plus.bedGraph chromInfo.txt $name2.plus.bw
	bedGraphToBigWig $name2.minus.bedGraph chromInfo.txt $name2.minus.bw
	
	rm *bed*	

	done
cd ..

###Generando BigWigs de alineamientos totales###
cd TOTAL

	for i in $(ls $name.TOTAL*.bam)

	do

	name3=${i%.bam}
	echo $name3

	bamToBed -i $name3.bam -split > $name3.bed

	awk '{if($6=="+") print}' $name3.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n > $name3.plus.bedGraph
	awk '{if($6=="-") print}' $name3.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > $name3.minus.bedGraph

	bedGraphToBigWig $name3.plus.bedGraph chromInfo.txt $name3.plus.bw
	bedGraphToBigWig $name3.minus.bedGraph chromInfo.txt $name3.minus.bw
	
	rm *bed*
	rm temporal*

	echo "Processing $name... OK" $(date)
	echo "Termino" $name | mail -s $name -a ../log.txt -c geparada88@gmail.com roberomunita@gmail.com

	done
cd ..

done


echo $(date)
echo "############################################################"
echo "#                                                          #"
echo "#                         DONE                             #"
echo "#                                                          #"
echo "############################################################"
