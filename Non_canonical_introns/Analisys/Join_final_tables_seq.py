import sys
import csv



def main(hg19, SJ_hg19, hg19_reads_seq_tags, GM12878, SJ_GM12878_paternal, SJ_GM12878_maternal, GM12878_reads_seq_tags_paternal, GM12878_reads_seq_tags_maternal, TOTAL_final_table):
	
	csv.field_size_limit(1000000000)
	reader1 = csv.reader(open(hg19), delimiter = ' ')
	reader2 = csv.reader(open(SJ_hg19), delimiter = ' ')
	reader3 = csv.reader(open(hg19_reads_seq_tags), delimiter = ' ')
	reader4 = csv.reader(open(GM12878), delimiter = ' ')
	reader5 = csv.reader(open(SJ_GM12878_paternal), delimiter = ' ')
	reader6 = csv.reader(open(SJ_GM12878_maternal), delimiter = ' ')
	reader7 = csv.reader(open(GM12878_reads_seq_tags_paternal), delimiter = ' ')
	reader8 = csv.reader(open(GM12878_reads_seq_tags_maternal), delimiter = ' ')
	reader9 = csv.reader(open(TOTAL_final_table), delimiter = ' ')
		
	hg19_intron_reads = {}
	GM12878_intron_reads = {}
	reads_seq = {}
	
	for row in reader1:
		intron = row[0]
		read = row[10].split(",")[:5]
		hg19_intron_reads[intron] = read
	
	for row in reader2:
		read = row[0]
		dn = row[7]
		seq = row[14]
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			reads_seq[read] = seq

	for row in reader3:
		read = row[0]
		seq = row[1]
		reads_seq[read] = seq
		
	for row in reader4:
		intron = row[0]
		read = row[10].split(",")[:5]
		GM12878_intron_reads[intron] = read			

	for row in reader5:
		read = row[0]
		dn = row[7]
		seq = row[14]
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			reads_seq[read] = seq

	for row in reader6:
		read = row[0]
		dn = row[7]
		seq = row[14]
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			reads_seq[read] = seq

	for row in reader7:
		read = row[0]
		seq = row[1]
		reads_seq[read] = seq
	
	for row in reader8:
		read = row[0]
		seq = row[1]
		reads_seq[read] = seq
	
	for row in reader9:
		intron = row[0]
		dn = row[6]
		hg19 = int(row[9])
		GM12878 = int(row[10])
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			
			seqs_hg19 = []
			seqs_GM12878 = []
			
			try:
				reads_hg19 = hg19_intron_reads[intron]
				for read in reads_hg19:
					seq = reads_seq[read]
					seqs_hg19.append(seq)
					
			except KeyError:
				pass
			
			try:
				reads_GM12878 = GM12878_intron_reads[intron]
				for read in reads_GM12878:
					seq = reads_seq[read]
					seqs_GM12878.append(seq)
					
			except KeyError:
				pass

			if seqs_hg19 == []:
				seqs_hg19 = "0"
			
			if seqs_GM12878 == []:
				seqs_GM12878 = "0"
			
			print " ".join(row), ",".join(seqs_hg19), ",".join(seqs_GM12878)
			
		
#python ~/my_src/Analisys/Join_final_tables_seq.py ../hg19/ALL/introns.final_table.hg19.fixed.tags ../hg19/ALL/SJ.introns.blat1.TOTAL tags/hg19/TOTAL.tags.filter.final ../GM12878/NA12878_Joel_Rozowsky/STRANDED/TOTAL/introns.final_table.hg19.fixed.tags ../GM12878/NA12878_Joel_Rozowsky/STRANDED/TOTAL_paternal/ALL/SJ.introns.blat1.TOTAL ../GM12878/NA12878_Joel_Rozowsky/STRANDED/TOTAL_maternal/ALL/SJ.introns.blat1.TOTAL tags/GM12878/paternal/TOTAL.tags.filter.final tags/GM12878/maternal/TOTAL.tags.filter.final TOTAL_introns.final_table.tags
	




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
