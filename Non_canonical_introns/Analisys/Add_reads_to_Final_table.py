import sys
import csv
from collections import defaultdict

def main(non_canonical_final_table, final_table, SJ_introns, reads_tags):
	csv.field_size_limit(1000000000)
	
	reader1 = csv.reader(open(non_canonical_final_table), delimiter = ' ')
	reader2 = csv.reader(open(final_table), delimiter = ' ')
	reader3 = csv.reader(open(SJ_introns), delimiter = ' ')
	reader4 = csv.reader(open(reads_tags), delimiter = ' ')

	
	Total_intron_info = {}
	non_can_reads_seq = {}
	introns_reads_tags = defaultdict(set)
	introns_seq_tags = defaultdict(set)
	
	data_set_introns = set([])
	
	for row in reader1:
		
		intron = row[0]
		chr = intron.split(":")[0]
		strand = ""
		if "+" in intron:
			strand = "+"
		elif "-" in intron:
			strand = "-"
		istart = intron.split(":")[1].split(strand)[0]
		iend = intron.split(":")[1].split(strand)[1]
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_score = row[8]		

		Total_intron_info[intron] = (intron, chr, strand, istart, iend, ilength, dn, dn_type, dn_score)
	
	for row in reader3:
		read = row[0]
		chr = row[1]
		istart = int(row[2])
		iend = int(row[3])
		strand = row[4]
		ilen = int(row[5])
		intron = row[6]
		dn = row[7]
		start = int(row[8])
		cigar = row[9]
		e5s = int(row[10])
		e5e = int(row[11])
		e3s = int(row[12])
		e3e = int(row[13])
		seq = row[14]
		end = e3e
		
		if dn!="GTAG" and dn!="GCAG" and dn!="ATAC":
			non_can_reads_seq[read] = seq
	
	for row in reader4:
		read = row[0]
		seq = row[1]
		qual = row[2]
		intron = row[3]
		up_anchor = int(row[4])
		down_anchor = int(row[5])
		
		introns_reads_tags[intron].add(read)
		introns_seq_tags[intron].add(seq)

	for row in reader2:
		intron = row[0]
		coverage = row[1]
		chr = row[2]
		strand = row[3]
		istart = row[4]
		iend = row[5]
		ilength = row[6]
		dn = row[7]
		dn_type = row[8]
		dn_score = row[9]
		reads = set(row[10].split(","))
		
		data_set_introns.add(intron)
		
		if dn=="GTAG" or dn=="GCAG" or dn=="ATAC":
			print " ".join(row)

		
		elif introns_reads_tags.has_key(intron):       #Esto es para los que ya estaban en la tabla final, pero se aumento el numero de reads encontrados
			seq_coverage = set([])
			for read in reads:
				seq = non_can_reads_seq[read]
				seq_coverage.add(seq)
			
			reads_tags = introns_reads_tags[intron]
			seq_tags = introns_seq_tags[intron]
			
			reads = reads | reads_tags
			seq_coverage = seq_coverage | seq_tags
			
			coverage = len(seq_coverage)
			
			print intron, coverage, chr, strand, istart, iend, ilength, dn, dn_type, dn_score, ",".join(reads)		 
		
		else:
			print " ".join(row)
	
	
	for row in introns_reads_tags.items():
		intron = row[0]
		reads = row[1]
		seq_coverage = introns_seq_tags[intron]
		
		
		if intron in data_set_introns:
			pass
		
		else:
			
			coverage = len(seq_coverage)

			chr = Total_intron_info[intron][1]
			strand = Total_intron_info[intron][2]
			istart = Total_intron_info[intron][3]
			iend = Total_intron_info[intron][4]
			ilength = Total_intron_info[intron][5]
			dn = Total_intron_info[intron][6]
			dn_type = Total_intron_info[intron][7]
			dn_score = Total_intron_info[intron][8]
			
			print intron, coverage, chr, strand, istart, iend, ilength, dn, dn_type, dn_score, ",".join(reads)

				

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
