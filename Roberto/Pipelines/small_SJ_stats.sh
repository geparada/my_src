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
