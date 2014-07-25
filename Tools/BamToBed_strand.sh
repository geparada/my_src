for i in $(ls *.bam)

do
	name=${i%.bam}
	echo $name

	bamToBed -i $name.bam -split > $name.bed

	awk '$6=="+" && NF==6' $name.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n > $name.plus.bedGraph &&
	awk '$6=="-" && NF==6' $name.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n > $name.minus.bedGraph &&

	bedGraphToBigWig $name.plus.bedGraph chromInfo.txt $name.plus.bw 
	bedGraphToBigWig $name.minus.bedGraph chromInfo.txt $name.minus.bw 

done
