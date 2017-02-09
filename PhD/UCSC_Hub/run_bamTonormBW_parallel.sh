#!/bin/bash

# samtools merge sRNA_SX1316.bam  /nfs/users/nfs_t/tdd/lustre/users/tdd/projects/emb4/out/sRNA_SX1316.?_seqfeatfilter_WS220-sensor.mis0.map5000.22G.bam 
# samtools merge sRNA_SX2000.bam  /nfs/users/nfs_t/tdd/lustre/users/tdd/projects/emb4/out/sRNA_SX2000.?_seqfeatfilter_WS220-sensor.mis0.map5000.22G.bam
# samtools merge sRNA_SX2929.bam  /nfs/users/nfs_t/tdd/lustre/users/tdd/projects/emb4/out/sRNA_SX2929.?_seqfeatfilter_WS220-sensor.mis0.map5000.22G.bam
# samtools merge sRNA_SX2930.bam  /nfs/users/nfs_t/tdd/lustre/users/tdd/projects/emb4/out/sRNA_SX2930.?_seqfeatfilter_WS220-sensor.mis0.map5000.22G.bam

samtools merge sRNA_SX2929-30.bam  /nfs/users/nfs_t/tdd/lustre/users/tdd/projects/emb4/out/sRNA_SX29??.?_seqfeatfilter_WS220-sensor.mis0.map5000.22G.bam


Genome="ce10"
fetchChromSizes $Genome > $Genome.chromsizes


for i in $(ls *.bam)


	do

		name=$(basename $i .bam)

		sed "s/NAME/$name/g" ~/my_src/PhD/UCSC_Hub/bamTonormBW.sh > bamTonormBW.$name.sh

		chmod +x bamTonormBW.$name.sh

<<<<<<< HEAD
		bash ~/submit.job -s ./bamTonormBW.$name.sh -m 80000 -q normal -n bw.$name.norm
=======
		bash ~/submit.job -s ./bamTonormBW.$name.sh -m 50000 -q normal -n bw.$name.norm
>>>>>>> 4c02a22aa98a8b48eb1592cab5f99c689f202615

	done
