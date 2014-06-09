import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


TAG_seq = {}
fasta = open("reads_tags.filter.fa", "w")

	
def TAGstabulator(tags_fasta):
	f = open(tags_fasta)
	for chrfa in SeqIO.parse(f, "fasta"):
		TAG_seq[str(chrfa.id)] = str(chrfa.seq)	
	f.close()


def main(tags_sam):
	reader1 = csv.reader(open(tags_sam), delimiter = '\t')
	
	row_1 = ["", "", ""]
	
	for row in reader1:
		if row_1[2] == "":
			print "\t".join(row)
		
		intron_1 = row_1[2].split("|")[0]
		read_1 = row_1[0]
		
		intron =  row[2].split("|")[0]
		read = row[0]
		

				
		
		
		if row[0][0]!='@' and (row[1]=='0' or row[1]=='16'):
			
			flag = int(row[1])
			tag_name = row[2]
			start = int(row[3])           #Sam es 1 referenciado 
			cigar = row[5]
			seq = row[9]
			qual = row[10]
			intron_tag = tag_name.split("|")[0]		
		
			mismatches = int(row[13].strip("NM:i"))
			matches = int(cigar.strip('M'))
			fist_block_tag = int(tag_name.split("|")[2].split("_")[0])
			anchor_up =  fist_block_tag - start
			anchor_down =  start + matches - fist_block_tag



			if read_1 != read:

				if mismatches <= 2 and anchor_up >= 8 and anchor_down >= 8:
					read_window = seq[fist_block_tag-start-8+1:fist_block_tag-start+8+1]
					seq_tag = TAG_seq[tag_name]
					tag_window = seq_tag[fist_block_tag-8:fist_block_tag+8]
					if read_window == tag_window:
						print "\t".join(row)
		
			elif intron_1 != intron:
			
				if mismatches <= 2 and anchor_up >= 8 and anchor_down >= 8:
					read_window = seq[fist_block_tag-start-8+1:fist_block_tag-start+8+1]
					seq_tag = TAG_seq[tag_name]
					tag_window = seq_tag[fist_block_tag-8:fist_block_tag+8]
					if read_window == tag_window:		
						print "\t".join(row)
		

		elif row[1]!='4':
			print "\t".join(row)


#		except IndexError:
#			print "\t".join(row)

#		except ValueError:
#			print "\t".join(row)		

#		row_1 = row
		


if __name__ == '__main__':
	TAGstabulator(sys.argv[1])
	main(sys.argv[2])
