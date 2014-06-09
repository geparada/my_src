import sys
from collections import defaultdict


def main (processed_list, trackDB_constant):
	""" Genera tracks de RNAs pequenos """
	
	for row in open(trackDB_constant, 'r'):
		print row.strip("\n")
	print "\n"
	
		
	cells_tracks = defaultdict(list)



	for row in open(processed_list, 'r'):
		file_base_name = row.strip(".fastq.clip.trim.gz\n")
		name = row.replace("wgEncodeCshlShortRnaSeq","").replace(".fastq.clip.trim.gz\n", "")
		
		cell_compartiment = ""
		if "Cell" in name:
			cell_compartiment = "Cell"

		if "Cytosol" in name:
			cell_compartiment = "Cytosol"
			
		if "Nucleus" in name:
			cell_compartiment = "Nucleus"
		
		cell_name = name.split(cell_compartiment)[0]
		
		treatement = "Total"
		if "ShorttotalCiptap" in name:
			treatement = "Cip Tap"
		if "ShorttotalTap" in name:
			treatement = "Tap"
		
		Rep = "Rep " + name.split("Rep")[-1]
		
		cells_tracks[cell_name].append((cell_compartiment, treatement, Rep, file_base_name))
		
	for key in cells_tracks.items():
		cell_name = key[0]

		BAM_SJ_ALL_ANTISENSE = []	
		BAM_SJ_UNIQUE_ANTISENSE = []
		BAM_SJ_ALL_SENSE = []		
		BAM_SJ_UNIQUE_SENSE = []
		
		BAM_SJ_ALL_ANTISENSE_FINAL_FILTERS = []	
		BAM_SJ_UNIQUE_ANTISENSE_FINAL_FILTERS = []
		BAM_SJ_ALL_SENSE_FINAL_FILTERS = []		
		BAM_SJ_UNIQUE_SENSE_FINAL_FILTERS = []
		
		BAM_polyTs = []
		BAM_polyTs_uniq = []		

		BW_SJ_ALL_PLUS= []	
		BW_SJ_ALL_MINUS = []
		BW_TOTAL_ALL_PLUS = []	
		BW_TOTAL_ALL_MINUS = []			
	
		BW_SJ_UNIQUE_PLUS = []
		BW_SJ_UNIQUE_MINUS= []		
		BW_TOTAL_UNIQUE_PLUS = []
		BW_TOTAL_UNIQUE_MINUS = []
		
		
		BAM_compositeTrack_name = cell_name + "_BAMs"
			
		
		for cell_info in key[1]:

			cell_compartiment, treatement, Rep, file_base_name = cell_info	

			BAM_polyTs.append(["track " + file_base_name + ".polyT.final", "parent BAM_polyTs"+ cell_name, "bigDataUrl " + file_base_name + ".polyT.final.bam", "shortLabel polyT Reads " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments", "longLabel PolyT Reads Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep, "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])
			BAM_polyTs_uniq.append(["track " + file_base_name + ".polyT.final.uniq", "parent BAM_polyTs_uniq"+ cell_name, "bigDataUrl " + file_base_name + ".polyT.final.uniq.bam", "shortLabel polyT Reads " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Uniq", "longLabel PolyT Reads Unique Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep, "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])
						
						
			BAM_SJ_ALL_ANTISENSE.append(["track " + file_base_name + ".antisense_SJ", "parent  BAM_SJ_ALL_ANTISENSE"+ cell_name, "bigDataUrl " + file_base_name + ".antisense_SJ.bam", "shortLabel Antisense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments", "longLabel Antisense Splice Junction Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep, "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])
			BAM_SJ_UNIQUE_ANTISENSE.append(["track " + file_base_name + ".antisense_SJ_uniq", "parent BAM_SJ_UNIQUE_ANTISENSE"+ cell_name, "bigDataUrl " + file_base_name + ".antisense_SJ_uniq.bam", "shortLabel Antisense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Uniq", "longLabel Antisense Splice Junction Unique Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep, "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])
			BAM_SJ_ALL_SENSE.append(["track " + file_base_name + ".sense_SJ", "parent BAM_SJ_ALL_SENSE"+ cell_name, "bigDataUrl " + file_base_name + ".sense_SJ.bam", "shortLabel Sense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments", "longLabel Sense Splice Junction Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep, "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])			
			BAM_SJ_UNIQUE_SENSE.append(["track " + file_base_name + ".sense_SJ_uniq", "parent BAM_SJ_UNIQUE_SENSE"+ cell_name, "bigDataUrl " + file_base_name + ".sense_SJ_uniq.bam", "shortLabel Sense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Uniq", "longLabel Sense Splice Junction Unique Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep, "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])

			BAM_SJ_ALL_ANTISENSE_FINAL_FILTERS.append(["track " + file_base_name + ".antisense_SJ.final_filters", "parent  BAM_SJ_ALL_ANTISENSE_FINAL_FILTERS"+ cell_name, "bigDataUrl " + file_base_name + ".antisense_SJ.final_filters.bam", "shortLabel Antisense SJ Final Filters " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments", "longLabel Antisense Splice Junction Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - Final filters"  , "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])
			BAM_SJ_UNIQUE_ANTISENSE_FINAL_FILTERS.append(["track " + file_base_name + ".antisense_SJ_uniq.final_filters", "parent BAM_SJ_UNIQUE_ANTISENSE_FINAL_FILTERS"+ cell_name, "bigDataUrl " + file_base_name + ".antisense_SJ_uniq.final_filters.bam", "shortLabel Antisense SJ Final Filters " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Uniq", "longLabel Antisense Splice Junction Unique Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - Final filters", "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])
			BAM_SJ_ALL_SENSE_FINAL_FILTERS.append(["track " + file_base_name + ".sense_SJ.final_filters", "parent BAM_SJ_ALL_SENSE_FINAL_FILTERS"+ cell_name, "bigDataUrl " + file_base_name + ".sense_SJ.final_filters.bam", "shortLabel Sense SJ Final Filters " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments", "longLabel Sense Splice Junction Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - Final filters", "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])			
			BAM_SJ_UNIQUE_SENSE_FINAL_FILTERS.append(["track " + file_base_name + ".sense_SJ_uniq.final_filters", "parent BAM_SJ_UNIQUE_SENSE_FINAL_FILTERS"+ cell_name, "bigDataUrl " + file_base_name + ".sense_SJ_uniq.final_filters.bam", "shortLabel Sense SJ Final Filters " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Uniq", "longLabel Sense Splice Junction Unique Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - Final filters", "type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])


			BW_SJ_ALL_PLUS.append(["track " + file_base_name + ".SJ.all.plus", "parent BW_SJ_ALL_PLUS"+ cell_name, "bigDataUrl " + file_base_name + ".SJ.all.plus.bw", "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Plus", "longLabel Coverage Splice Junction of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand", "type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])
			BW_SJ_ALL_MINUS.append(["track " + file_base_name + ".SJ.all.minus", "parent BW_SJ_ALL_MINUS"+ cell_name, "bigDataUrl " + file_base_name + ".SJ.all.minus.bw", "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Minus", "longLabel Coverage Splice Junction of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand", "type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])
			BW_TOTAL_ALL_PLUS.append(["track " + file_base_name + ".TOTAL.all.plus", "parent BW_TOTAL_ALL_PLUS"+ cell_name, "bigDataUrl " + file_base_name + ".TOTAL.all.plus.bw", "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Plus", "longLabel Coverage of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand", "type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])
			BW_TOTAL_ALL_MINUS.append(["track " + file_base_name + ".TOTAL.all.minus", "parent  BW_TOTAL_ALL_MINUS"+ cell_name, "bigDataUrl " + file_base_name + ".TOTAL.all.minus.bw", "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Minus", "longLabel Coverage of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand", "type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])		

			BW_SJ_UNIQUE_PLUS.append(["track " + file_base_name + ".SJ.uniq.plus", "parent BW_SJ_UNIQUE_PLUS"+ cell_name, "bigDataUrl " + file_base_name + ".SJ.uniq.plus.bw", "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Plus", "longLabel Coverage Splice Junction of Unique Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand", "type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			BW_SJ_UNIQUE_MINUS.append(["track " + file_base_name + ".SJ.uniq.minus", "parent BW_SJ_UNIQUE_MINUS"+ cell_name, "bigDataUrl " + file_base_name + ".SJ.uniq.minus.bw", "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Minus", "longLabel Coverage Splice Junction of Unique Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand", "type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			BW_TOTAL_UNIQUE_PLUS.append(["track " + file_base_name + ".TOTAL.uniq.plus", "parent BW_TOTAL_UNIQUE_PLUS"+ cell_name, "bigDataUrl " + file_base_name + ".TOTAL.uniq.plus.bw", "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Plus", "longLabel Coverage of Unique Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand", "type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])
			BW_TOTAL_UNIQUE_MINUS.append(["track " + file_base_name + ".TOTAL.uniq.minus", "parent BW_TOTAL_UNIQUE_MINUS"+ cell_name, "bigDataUrl " + file_base_name + ".TOTAL.uniq.minus.bw", "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Minus", "longLabel Coverage of Unique Alignments" + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand", "type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			


				########## BAMs ##########

		BAMs_polyT_compositeTrack_name = cell_name + "polyT_BAMs"	

		print "\t" + "track " + BAMs_polyT_compositeTrack_name
		print "\t" + "shortLabel " + cell_name + " polyT Reads Alignments"
		print "\t" + "longLabel Alignments of polyT Reads " + cell_name  				
		print "\t" + "superTrack on"
		print "\n"	

		print "\t\t" + "track BAM_polyTs" + cell_name
		print "\t\t" + "parent " + BAMs_polyT_compositeTrack_name
		print "\t\t" + "shortLabel " + cell_name + " PolyT Reads"
		print "\t\t" + "longLabel Alignments of Reads With PolyT - " + cell_name  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_polyTs:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"


		print "\t\t" + "track BAM_polyTs_uniq" + cell_name
		print "\t\t" + "parent " + BAMs_polyT_compositeTrack_name
		print "\t\t" + "shortLabel " + cell_name + " PolyT Unique Reads"
		print "\t\t" + "longLabel Unique Alignments of Reads With PolyT - " + cell_name  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_polyTs_uniq:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"


####### SJ
			
		BAMs_SJ_compositeTrack_name = cell_name + "SJ_BAMs"			

		print "\t" + "track " + BAMs_SJ_compositeTrack_name
		print "\t" + "shortLabel " + cell_name + " Spliced Alignments"
		print "\t" + "longLabel Alignments of Spliced smallRNAs - " + cell_name  				
		print "\t" + "superTrack on"
		print "\n"	


		
		print "\t\t" + "track BAM_SJ_ALL_ANTISENSE" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Antisense Alignments"
		print "\t\t" + "longLabel Antisense Alignments of Spliced smallRNAs - " + cell_name  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_ALL_ANTISENSE:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"
				
		print "\t\t" + "track BAM_SJ_ALL_SENSE" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Sense Alignments"
		print "\t\t" + "longLabel Alignments of Sense Spliced  smallRNAs - " + cell_name  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_ALL_SENSE:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BAM_SJ_UNIQUE_ANTISENSE" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Antisense Unique Alignments"
		print "\t\t" + "longLabel Unique Antisense Alignments of Spliced smallRNAs - " + cell_name  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_UNIQUE_ANTISENSE:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"


		print "\t\t" + "track BAM_SJ_UNIQUE_SENSE" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Sense Unique Alignments"
		print "\t\t" + "longLabel Unique Alignments of Sense Spliced  smallRNAs - " + cell_name  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_UNIQUE_SENSE:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"





###########
		
		print "\t\t" + "track BAM_SJ_ALL_ANTISENSE_FINAL_FILTERS" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Antisense Alignments Final Filters"
		print "\t\t" + "longLabel Antisense Alignments of Spliced smallRNAs - " + cell_name + " - Final Filters"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_ALL_ANTISENSE_FINAL_FILTERS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"
				
		print "\t\t" + "track BAM_SJ_ALL_SENSE_FINAL_FILTERS" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Sense Alignments Final Filters"
		print "\t\t" + "longLabel Alignments of Sense Spliced  smallRNAs - " + cell_name + " - Final Filters"   				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_ALL_SENSE_FINAL_FILTERS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BAM_SJ_UNIQUE_ANTISENSE_FINAL_FILTERS" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Antisense Unique Alignments Final Filters"
		print "\t\t" + "longLabel Unique Antisense Alignments of Spliced smallRNAs - " + cell_name + " - Final Filters"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_UNIQUE_ANTISENSE_FINAL_FILTERS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"


		print "\t\t" + "track BAM_SJ_UNIQUE_SENSE_FINAL_FILTERS" + cell_name
		print "\t\t" + "parent " + BAMs_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Sense Unique Alignments Final Filters"
		print "\t\t" + "longLabel Unique Alignments of Sense Spliced  smallRNAs - " + cell_name + " - Final Filters"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BAM_SJ_UNIQUE_SENSE_FINAL_FILTERS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"



###########

		########## bws ##########

		bws_SJ_compositeTrack_name = cell_name + "SJ_bws"

		print "\t" + "track " + bws_SJ_compositeTrack_name
		print "\t" + "shortLabel " + cell_name + " Coverage"
		print "\t" + "longLabel Coverage of smallRNAs - " + cell_name   				
		print "\t" + "superTrack on"
		print "\n"	

		print "\t\t" + "track BW_SJ_ALL_PLUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name +  " SJ Coverage Plus"
		print "\t\t" + "longLabel Coverage of Spliced smallRNAs - " + cell_name + " - Plus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_SJ_ALL_PLUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BW_SJ_ALL_MINUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Coverage Minus"
		print "\t\t" + "longLabel Coverage of Spliced smallRNAs - " + cell_name + " - Minus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_SJ_ALL_MINUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BW_TOTAL_ALL_PLUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " Total Coverage Plus"
		print "\t\t" + "longLabel Coverage of smallRNAs - " + cell_name + " - Plus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_TOTAL_ALL_PLUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BW_TOTAL_ALL_MINUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " Total Coverage Minus"
		print "\t\t" + "longLabel Coverage of smallRNAs - " + cell_name + " - Minus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_TOTAL_ALL_MINUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"


		print "\t\t" + "track BW_SJ_UNIQUE_PLUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Unique Coverage Plus"
		print "\t\t" + "longLabel Unique Coverage of Spliced smallRNAs - " + cell_name + " - Plus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_SJ_UNIQUE_PLUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BW_SJ_UNIQUE_MINUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " SJ Unique Coverage Minus"
		print "\t\t" + "longLabel Unique Coverage of Spliced smallRNAs - " + cell_name + " - Minus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_SJ_UNIQUE_MINUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BW_TOTAL_UNIQUE_PLUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " Total Unique Coverage Plus"
		print "\t\t" + "longLabel Unique Coverage of smallRNAs - " + cell_name + " - Plus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_TOTAL_UNIQUE_PLUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"

		print "\t\t" + "track BW_TOTAL_UNIQUE_MINUS" + cell_name
		print "\t\t" + "parent " + bws_SJ_compositeTrack_name
		print "\t\t" + "shortLabel Spliced " + cell_name + " Total Unique Coverage Minus"
		print "\t\t" + "longLabel Unique Coverage of smallRNAs - " + cell_name + " - Minus"  				
		print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
		print "\n"	

		for data_file in BW_TOTAL_UNIQUE_MINUS:
			for track_line in data_file:
				print "\t\t\t" + track_line
			print "\n"




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
