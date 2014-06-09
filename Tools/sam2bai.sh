for i in $(ls *.sam)

do

name=${i%.sam}
samtools view -Sb $name.sam -o $name.bam && samtools sort $name.bam $name.sort.bam $name.sort.bam.bai && rm $name.sam $name.bam

done
