#Para que este script funcione, los reads alineados deben esta en una carpeta del mismo nombre que su FASTQ
# Ademas se deve correr desde la carpeta que contiene todas las carpetas de los reads alineados

#Completa las variables

Genome="/home/geparada/db/genome/hg19.fa"
READS_NAMES="_clip*"               #Exprecion Regular que contiene los distindos grupos de reads
FORWARD_PAIR="Rd1"                 #Poner cual de los dos pares se secuencio forward Rd1 o Rd2
ANCHOR="8"
COVERAGE="3"
MIN_INTRON_LENGTH="6"          #Este es el rango en que MapSplice 2.0 toma los gaps como intrones...
MAX_INTRON_LENGTH="200000"     

#CHAIN_FILE="/home/geparada/db/genome/GM12878/Joel_Rozowsky/paternal.chain.swap"



#Pipeline:

echo "Extracting introns from SAM:"

mkdir ALL && mkdir Repbase

for i in $(ls -d $READS_NAMES)

	do
	 
	cd $i &&
#	samtools view alignments.bam -H > alignments.sam.only-uniq
#	samtools view alignments.bam | grep -e IH:i:1 -w >> alignments.sam.only-uniq && 
#	python ~/my_src/Tools/Get_introns_from_sam.py $Genome alignments.sam.only-uniq $FORWARD_PAIR $MIN_INTRON_LENGTH $MAX_INTRON_LENGTH $ANCHOR > SJ.introns &&               #viene incluido el filtro por anchor
#	samtools view -S -b -h alignments.sam.only-uniq -o alignments.bam.only-uniq &&
#	rm alignments.sam.only-uniq &&
	cd ..

	done

echo "NON-CANONICAL BLAT FILTER:"

for i in $(ls -d $READS_NAMES)

	do
	 
	echo Proccesing $i ...  && 
	cd $i &&
#	python ~/my_src/Tools/FASTA_from_non_canonical_spliced_reads.py SJ.introns > non_canonical.fasta &&         
#	dust non_canonical.fasta | grep NNNNN -B 1 | grep ^\>  >> ../Repbase/low_complexity_IDs.dust  &&
#	blat $Genome non_canonical.fasta non_canonical.psl -noHead &&         # Los reads que tienen intrones no cacnonicos se re-mapean con blat
#	python ~/my_src/Filters/Error_classifier.py $Genome SJ.introns non_canonical.psl > SJ.introns.blat1_error &&   #Se identifican los reads que hacen multimapping, que alinean en otro locus segunn BLAT o que al correr el alineamiento por DRs quedan como intrones canonicos
#	grep -v -e BAD_ALIGNMENT SJ.introns.blat1_error > SJ.introns.blat1  && 
	cd ..

	done

echo "MERGING SPLICED READS:"

#cat $READS_NAMES/non_canonical.fasta > Repbase/non_canonical.TOTAL.fasta
#cat $READS_NAMES/SJ.introns.blat1 > ALL/SJ.introns.blat1.TOTAL
#cat $READS_NAMES/non_canonical.psl > ALL/non_canonical_blat_genome


echo "SEARCHING SPLICED READS ON REPBASE:"

cd Repbase && 
#blat ~/RepeatMasker/Libraries/Homo_sapiens_all_Repbase non_canonical.TOTAL.fasta reads.Repbase.psl -fastMap -noHead && 
cd ..

echo "FINAL FILTERS:"

#python ~/my_src/Filters/Filtro_anchor_repbase_dust_coverage.py $Genome ALL/SJ.introns.blat1.TOTAL Repbase/reads.Repbase.psl Repbase/low_complexity_IDs.dust $COVERAGE > ALL/introns.pre_final_table &&
##liftOver ALL/introns.pre_final_table $CHAIN_FILE ALL/introns.pre_final_table.hg18 ALL/introns.pre_final_table.hg18.unmapped  &&
##liftOver ALL/introns.pre_final_table.hg18 ~/db/genome/liftover_fies/hg18Tohg19.over.chain ALL/introns.pre_final_table.hg19 ALL/introns.pre_final_table.hg19.unmapped &&
#sed 's/ /\t/g' ALL/introns.pre_final_table > ALL/introns.pre_final_table.hg19
##liftOver ALL/non_canonical_blat_genome $CHAIN_FILE ALL/non_canonical_blat_genome.hg18 ALL/non_canonical_blat_genome.hg18.unmapped -pslT &&
##liftOver ALL/non_canonical_blat_genome.hg18 ~/db/genome/liftover_fies/hg18Tohg19.over.chain ALL/non_canonical_blat_genome.hg19 ALL/non_canonical_blat_genome.hg19.unmapped -pslT &&
#mv ALL/non_canonical_blat_genome ALL/non_canonical_blat_genome.hg19
#python ~/my_src/Tools/non_canonical_SJTags.py ~/db/genome/hg19.fa ~/db/transcriptome/hg19/Gene_models/gencode/v11/gencode.v11.pc_transcripts.fa ~/db/transcriptome/hg19/Gene_models/gencode/v11/gencode.v11.annotation.bed12 ALL/introns.pre_final_table.hg19 > ALL/non_canonical_tags.fasta &&
#blat ALL/non_canonical_tags.fasta Repbase/non_canonical.TOTAL.fasta ALL/non_canonical_tags.psl -noHead -fastMap &&
python ~/my_src/Filters/Final_BLAT_filter.py ~/db/genome/hg19.fa ALL/introns.pre_final_table.hg19 ALL/non_canonical_blat_genome.hg19 ALL/non_canonical_tags.psl > ALL/introns.final_table.hg19 &&
python ~/my_src/Tools/Define_dn_PWM.py ~/db/genome/hg19.fa ALL/introns.final_table.hg19 /home/geparada/db/PWM/hg19_GT_AG_U2_5.good.matrix /home/geparada/db/PWM/hg19_GT_AG_U2_3.good.matrix /home/geparada/db/PWM/hg19_GC_AG_U2_5.good.matrix /home/geparada/db/PWM/hg19_GC_AG_U2_3.good.matrix /home/geparada/db/PWM/hg19_GT_AG_U12_5.good.matrix /home/geparada/db/PWM/hg19_GT_AG_U12_3.good.matrix /home/geparada/db/PWM/hg19_AT_AC_U12_5.good.matrix /home/geparada/db/PWM/hg19_AT_AC_U12_3.good.matrix > ALL/introns.final_table.hg19.fixed &&
awk '{if ($7>=40) print $0}' ALL/introns.final_table.hg19.fixed > ALL/introns.final_table.hg19.fixed.longer40pb &&
python ~/my_src/Tools/Dinucleotide_Pecentage.py ALL/introns.final_table.hg19.fixed.longer40pb > ALL/introns.final_table.hg19.fixed.longer40pb.dn.frec &&
#cat ALL/introns.final_table.hg19.fixed.longer40pb.dn.frec | mail -s "INTRONES DE Illumina BodyMap 2.0 PROCESADOS EXITOSAMENTE" geparada88@gmail.com -c robertomunita@gmail.com


#OJO que las pre tablas finales despues de hacer el liftOver quedan separadas por tabular, asi que al adaptar el codigo al ilumina BodyMap, hay que reemplazar los liftOvers por sed 's/ /\t/g'
