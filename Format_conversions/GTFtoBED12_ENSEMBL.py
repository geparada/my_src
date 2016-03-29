import sys
import csv
from collections import defaultdict



keys = []
exon_start = []
exon_size = []
transcript_start_end = []

transcript_id_info = []

def converter (GTF):
	#reader = csv.reader(open(GTF), delimiter = '\t')
	
	f = open(GTF)
	
	
	for line in f:
		row = line.replace(";","\t").replace(" ","\t").replace("\t\t","\t").split("\t")		
		try:

			#print row[2]
			chr = row[0]
			#print row
			group = row[1]
			blocktype = row[2]
			block_start = int(row[3]) - 1
			block_end = int(row[4])
			block_size = block_end - block_start
			strand = row[6]
			gene_id = row[8]
			transcript_id = row[13]
			gene_type = row[23]
			#print row
			gene_status = row[15]
			gene_name = row[19]
			transcript_type = row[23]
			transcript_status = row[21]
			transcript_name = row[25]
			#print row
			#level = int(row[25])
			level = 0
			
			#print row[:3]
			#print blocktype, transcript_id, gene_type, transcript_name, transcript_type

			if blocktype == 'exon':
				#print transcript_id, gene_type, transcript_name, transcript_type
				keys.append(transcript_id)
				transcript_id_info.append((transcript_id, [chr, group, blocktype, strand, gene_id, gene_type, gene_status, gene_name, transcript_type, transcript_status, transcript_name, level]))
						
				exon_start.append((transcript_id, block_start))
				exon_size.append((transcript_id, block_size))
				
			if blocktype == 'transcript':
				#print transcript_id, gene_type, transcript_name, transcript_type
				transcript_start_end.append((transcript_id, [block_start, block_end]))
				
			
		except IndexError:             
			pass
		
		except ValueError:
			pass

	dict_id_info = dict(transcript_id_info)
	dict_transcript_start_end = dict(transcript_start_end)

#	print dict_transcript_start_end

	dict_exon_start = defaultdict(list)
	for transcript_id, block_start in exon_start:
		dict_exon_start[transcript_id].append(block_start)

	dict_exon_size = defaultdict(list)
	for transcript_id, block_size in exon_size:
		dict_exon_size[transcript_id].append(block_size)

	for key in set(keys):		

		chr = dict_id_info[key][0]
		name = key
		strand = dict_id_info[key][3]
		blocknum = len(dict_exon_size[key])
		blocksizes = dict_exon_size[key]
		exon_starts = dict_exon_start[key]
		if strand == '-':
			blocksizes = dict_exon_size[key][::-1]
			exon_starts = dict_exon_start[key][::-1]
		start = dict_transcript_start_end[key][0]
		end = dict_transcript_start_end[key][1]
		tstarts = [ x-start for x in exon_starts]
				
		group = dict_id_info[key][1]
		blocktype = dict_id_info[key][2]
		gene_id = dict_id_info[key][4]
		gene_type = dict_id_info[key][5]
		gene_status = dict_id_info[key][6]
		gene_name = dict_id_info[key][7]
		transcript_type = dict_id_info[key][8]
		transcript_status = dict_id_info[key][9]
		transcript_name = dict_id_info[key][10]
		level = dict_id_info[key][11]

		###All-info###
		BED12 = [chr, str(start), str(end), name.strip('"'), gene_name.strip('"'), strand, transcript_name.strip('"'), transcript_type.strip('"'), str(level), str(blocknum), ",".join(map(str, blocksizes)) + "," , ",".join(map(str, tstarts)) + ",", group, gene_id.strip('"'), gene_status.strip('"'), transcript_status.strip('"'), gene_type.strip('"') ]

		###Standar###
		BED12 = [chr, str(start), str(end), name.strip('"'), "0", strand, str(start), str(end), "255,0,0", str(blocknum), ",".join(map(str, blocksizes)) + "," , ",".join(map(str, tstarts)) + ","] 

		print "\t".join(BED12)		

		#chr, start, end, name, gene_name, strand, transcript_name, transcript_type, level, blocknum, blocksizes, tstarts, group, gene_id, gene_status, transcript_status, gene_type

		


#		print key, dict_exon_size[key], dict_exon_start[key], dict_id_info[key] 
	

	f.close()








if __name__ == '__main__':
	converter (sys.argv[1])
		
		
