
for i in $files

	do
	 
	name=${i%.fastq.gz}
	mv $name.fastq.gz

	parallel --trc {.}.out: \
"fastx_clipper -a AAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -l 22 -i {} -o {.}.out"




#cat list_of_files | \
#parallel --trc {.}.out -S server1,server2,: \
#"my_script {} > {.}.out"
