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
		
		super_track_name = "smallRNAs_" + cell_name + "_total"
		
		print "track " + super_track_name
		print "superTrack on"
		print "shortLabel SmallRNAs " + cell_name
		print "longLabel SmallRNAs on " + cell_name + " cell line"
		print "\n"
		
		for cell_info in key[1]:

			cell_compartiment, treatement, Rep, file_base_name = cell_info
			
			print "\t### ---- Tracks " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " ---- ###"
			print "\n"			
			
		
			
			########## BAMs ##########
			
			BAM_compositeTrack_name = ("_").join([cell_name, cell_compartiment, treatement, Rep]).replace(" ", "_") + "_BAMs"			

			print "\t\t" + "track " + BAM_compositeTrack_name
			print "\t\t" + "parent " + super_track_name
			print "\t\t" + "shortLabel Spliced " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments"
			print "\t\t" + "longLabel Alignments of Spliced " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep   				
			print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
			print "\n"	

			print "\t\t\t" + "track " + file_base_name + ".antisense_SJ"
			print "\t\t\t" + "parent " + BAM_compositeTrack_name
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".antisense_SJ.bam"
			print "\t\t\t" + "shortLabel Antisense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments"			
			print "\t\t\t" + "longLabel Antisense Splice Junction Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep
			print "\t\t\t" + "\n\t\t\t".join(["type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])						
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".sense_SJ"
			print "\t\t\t" + "parent " + BAM_compositeTrack_name
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".sense_SJ.bam"
			print "\t\t\t" + "shortLabel Sense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments"			
			print "\t\t\t" + "longLabel Sense Splice Junction Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep
			print "\t\t\t" + "\n\t\t\t".join(["type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])						
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".antisense_SJ_uniq"
			print "\t\t\t" + "parent " + BAM_compositeTrack_name
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".antisense_SJ_uniq.bam"
			print "\t\t\t" + "shortLabel Antisense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Uniq"			
			print "\t\t\t" + "longLabel Antisense Splice Junction Unique Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep
			print "\t\t\t" + "\n\t\t\t".join(["type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])						
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".sense_SJ_uniq"
			print "\t\t\t" + "parent " + BAM_compositeTrack_name
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".sense_SJ_uniq.bam"
			print "\t\t\t" + "shortLabel Sense SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Uniq"			
			print "\t\t\t" + "longLabel Sense Splice Junction Unique Alignments of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep
			print "\t\t\t" + "\n\t\t\t".join(["type bam", "showNames on", "bamColorMode strand", "indelDoubleInsert on", "indelQueryInsert on", "maxWindowToDraw 1000000"])						
			print "\n"

			########## bws ##########

			bw_compositeTrack_name = ("_").join([cell_name, cell_compartiment, treatement, Rep]).replace(" ", "_") + "_bws"

			print "\t\t" + "track " + bw_compositeTrack_name
			print "\t\t" + "parent " + super_track_name
			print "\t\t" + "shortLabel " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Coverage"
			print "\t\t" + "longLabel Coverage of " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep   				
			print "\t\t" + "\n\t\t".join(["compositeTrack on", "type bam 0 1.0", "viewLimits 0.0:0.2", "allButtonPair on", "dragAndDrop on", "centerLabelsDense on", "showSubtrackColorOnUi on"])
			print "\n"	

			print "\t\t\t" + "track " + file_base_name + ".SJ.all.plus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".SJ.all.plus.bw"
			print "\t\t\t" + "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Plus"						
			print "\t\t\t" + "longLabel Coverage Splice Junction of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".SJ.all.minus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".SJ.all.minus.bw"
			print "\t\t\t" + "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Minus"						
			print "\t\t\t" + "longLabel Coverage Splice Junction of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".TOTAL.all.plus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".TOTAL.all.plus.bw"
			print "\t\t\t" + "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Plus"						
			print "\t\t\t" + "longLabel Coverage of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".TOTAL.all.minus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".TOTAL.all.minus.bw"
			print "\t\t\t" + "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Alignments Minus"						
			print "\t\t\t" + "longLabel Coverage of Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".SJ.uniq.plus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".SJ.uniq.plus.bw"
			print "\t\t\t" + "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Plus"						
			print "\t\t\t" + "longLabel Coverage Splice Junction of Unique Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".SJ.uniq.minus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".SJ.uniq.minus.bw"
			print "\t\t\t" + "shortLabel Coverage SJ " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Minus"						
			print "\t\t\t" + "longLabel Coverage Splice Junction of Unique Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".TOTAL.uniq.plus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".TOTAL.uniq.plus.bw"
			print "\t\t\t" + "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Plus"						
			print "\t\t\t" + "longLabel Coverage of Unique Alignments " + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Plus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 0,0,150", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"

			print "\t\t\t" + "track " + file_base_name + ".TOTAL.uniq.minus"
			print "\t\t\t" + "parent " + bw_compositeTrack_name			
			print "\t\t\t" + "bigDataUrl " + file_base_name + ".TOTAL.uniq.minus.bw"
			print "\t\t\t" + "shortLabel Coverage " + (" ").join([cell_name, cell_compartiment, treatement, Rep]) + " Unique Alignments Minus"						
			print "\t\t\t" + "longLabel Coverage of Unique Alignments" + treatement + " smallRNAs - " + cell_name + "/" + cell_compartiment + " - " + Rep + " - " + "Minus Strand"
			print "\t\t\t" + "\n\t\t\t".join(["type bigWig", "color 150,0,0", "configurable on", "maxHeightPixels 30:30:30", "autoScale on"])			
			print "\n"





if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
