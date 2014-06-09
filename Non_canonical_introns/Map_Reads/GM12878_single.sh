
# Usage:
#     mapsplice [Inputs] [options]
#
# Inputs:
#    -1/--end_1                     <string>    [ end 1 reads or single end reads          ]
#    -2/--end_2                     <string>    [ end 2 reads                              ]
#    -c/--chromosome-dir            <string>    [ input chromosomes directory              ]
#    -x                             <string>    [ path and prefix of bowtie index          ]
# Options:
#    -p/--threads                   <int>       [ default: 1                               ]
#    -o/--output                    <string>    [ default: ./mapsplice_out                 ]
#    --bam                                      [ generate bam output                      ]
#    --keep-tmp                                 [ keep the intermediate files              ]
#    -l/--seglen                    <int>       [ default: 25                              ]
#    --min-map-len                  <int>       [ default: 0                               ]
#    --min-len                      <int>       [ default: 25                              ]
#    -m/--splice-mis                <int>       [ default: 1                               ]
#    --max-append-mis               <int>       [ default: 3                               ]
#    -k/max-hits                    <int>       [ default: 4                               ]
#    -i/--min-intron                <int>       [ default: 50                              ]
#    -I/--max-intron                <int>       [ default: 200,000                         ]
#    --qual-scale                   <string>    [ phred64(default) or phred33 or solexa64  ]
#    --del                          <int>       [ default: 6                               ]
#    --ins                          <int>       [ default: 6                               ]
#    --non-canonical                            [ output noncanonical junction             ]
#    --fusion | --fusion-non-canonical          [ output fusion junction                   ]
#    -h/--help                                  [ print the usage message                  ]
#    -v/--version                               [ print the version of MapSplice           ]


# Ejecutarse como bash GINGERAS_paired.sh
# En el READS_PATH deben haber dos carpetas gemelas llamadas Rd1 y Rd2  

Genome="/home/geparada/db/genome/GM12878/Joel_Rozowsky/paternal/"
READS_PATH="/home/geparada/db/transcriptome/GM12878/GINGERAS/single_end"            
READS_NAMES="wgEncodeCshlLongRnaSeqGm12878C*"                #Expresion regular que agrupe a los reads
INDEX="bowtie_index_GM12878_Joel_Rozowsky_paternal"
OUT_PATH="/media/HD2/Resultados/GM12878/NA12878_Joel_Rozowsky/GINGERAS/single_end/paternal/"
QUAL_SCALE="phred64"                                        #phred64 or phred33 or solexa64
   

cd ~/MapSplice_multi_threads_2.0_beta

for i in $(ls $READS_PATH/$READS_NAMES)

	do

	PRE_OUT_NAME="${i/$READS_PATH//}"
	OUT_NAME="${PRE_OUT_NAME////}"

	echo $OUT_NAME

	echo "python bin/mapsplice_multi_thread.py -1 $i -c $Genome -x $INDEX -o $OUT_PATH$OUT_NAME -p 8 -non-canonical --qual-scale $QUAL_SCALE"

	python bin/mapsplice_multi_thread.py -1 $i -c $Genome -x $INDEX -o $OUT_PATH$OUT_NAME -p 8 -non-canonical --qual-scale $QUAL_SCALE 

	done


Genome="/home/geparada/db/genome/GM12878/Joel_Rozowsky/maternal/"
READS_PATH="/home/geparada/db/transcriptome/GM12878/GINGERAS/single_end"            
READS_NAMES="wgEncodeCshlLongRnaSeqGm12878C*"                #Expresion regular que agrupe a los reads
INDEX="bowtie_index_GM12878_Joel_Rozowsky_maternal"
OUT_PATH="/media/HD2/Resultados/GM12878/NA12878_Joel_Rozowsky/GINGERAS/single_end/maternal/"
QUAL_SCALE="phred64"                                        #phred64 or phred33 or solexa64
   

cd ~/MapSplice_multi_threads_2.0_beta

for i in $(ls $READS_PATH/$READS_NAMES)

	do

	PRE_OUT_NAME="${i/$READS_PATH//}"
	OUT_NAME="${PRE_OUT_NAME////}"

	echo $OUT_NAME

	echo "python bin/mapsplice_multi_thread.py -1 $i -c $Genome -x $INDEX -o $OUT_PATH$OUT_NAME -p 8 -non-canonical --qual-scale $QUAL_SCALE"

	python bin/mapsplice_multi_thread.py -1 $i -c $Genome -x $INDEX -o $OUT_PATH$OUT_NAME -p 8 -non-canonical --qual-scale $QUAL_SCALE 

	done



echo 'MAPSPLICE 2.0 TERMINO DE MAPPEAR' | mail -s 'MapSplice' geparada88@gmail.com







