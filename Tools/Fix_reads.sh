fastq='10516_30DY0AAXX_c151_l6_r2.fastq'

#for i in $(ls 10516_30DY0AAXX_c151_l6_r2.fastq):

#do

awk 'BEGIN{OFS="\n"} { 
        a[NR % 4] = $0; 
        if(NR % 4 == 0 && length(a[2]) == length(a[0])){
            print a[1],a[2],a[3],a[0] 
        }
    }' $fastq   > $fastq.fixed

#done
