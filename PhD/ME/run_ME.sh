
#!/bin/bash

file_path="/lustre/scratch108/compgen/team218/gp7/Micro-exons/Mouse_ENCODE/fastq"
STRANDED="F"

TAGs="/lustre/scratch108/compgen/team218/gp7/Genome/mm10/mm10.ME_TAGs.fa"






# file_path="/lustre/scratch108/compgen/team218/gp7/Micro-exons/IBM2.0"
# STRANDED="T-F-DONE"

# TAGs="/lustre/scratch108/compgen/team218/gp7/Genome/hg19/hg19.ME_TAGs"




for i in $(ls $file_path/*.fastq.gz)

#for i in $(ls $file_path/ERR0308{56..71}.sra) #Mixed - standed

#for i in $(ls $file_path/ERR030{872..903}.sra) #Tissues - unstranded

 	do 

 		name=$(basename $i .fastq.gz)

 		echo $name


		# while [ $(bjobs | grep gp7 | wc -l ) -gt 1 ]
		# 	do
		# 		N=$(bjobs | grep gp7 | wc -l)
		# 		T=25

		# 		if [ $N -gt $T ]
		# 			then
		# 				echo OVER THRESHOLD \($T JOBS\)

		# 			else

		# 				break
		# 		fi

		# 		sleep 1
		# 	done


 		# sed "s/NAME/$name/g" /nfs/users/nfs_g/gp7/my_src/PhD/ME/Round1.sh | sed "s/STRANDED/$STRANDED/g" > $name.round1.sh
 		# chmod +x $name.round1.sh
 		# bash ~/submit.job -s ./$name.round1.sh -m 30000 -n R1.$name -q normal


 		sed "s/NAME/$name/g" /nfs/users/nfs_g/gp7/my_src/PhD/ME/Round2.sh | sed "s/STRANDED/$STRANDED/g" > $name.round2.sh
 		chmod +x $name.round2.sh
 		bash ~/submit.job -s ./$name.round2.sh -m 30000 -n R1.$name -q normal


	 done







