for i in $(ls *.bam)

do
	name=${i%.bam}
	echo $name

	bamToBed -i $name.bam -split > $name.bed
	awk 'NF==6' $name.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=chromInfo.txt stdin | sort -k1,1 -k2,2n > $name.bedGraph
	bedGraphToBigWig $name.bedGraph chromInfo.txt $name.bw &&
	rm $name.bed $name.bedGraph

done
