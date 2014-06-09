import sys
import csv
import pysam
from decimal import *
getcontext().prec = 20
csv.field_size_limit(1000000000)
import re

#Este programa es para asignarle un fold a los intron retention

#M 0 alignment match (can be a sequence match or mismatch)
#I 1 insertion to the reference
#D 2 deletion from the reference
#N 3 skipped region from the reference
#S 4 soft clipping (clipped sequences present in SEQ)
#H 5 hard clipping (clipped sequences NOT present in SEQ)
#P 6 padding (silent deletion from padded reference)
#= 7 sequence match
#X 8 sequence mismatch

def main(Final_table, paternal_pre_final_table, maternal_pre_final_table, bam, paternal_bam, maternal_bam):
	

	reader1 = csv.reader(open(Final_table), delimiter = '\t')
	reader2 = csv.reader(open(paternal_pre_final_table), delimiter = '\t')
	reader3 = csv.reader(open(maternal_pre_final_table), delimiter = '\t')
	hg19_bam = pysam.Samfile( bam, "rb" )
	paternal_bam = pysam.Samfile( paternal_bam, "rb" )
	maternal_bam = pysam.Samfile( maternal_bam, "rb" )
	
	hg19_to_paternal = {}
	hg19_to_maternal = {}
	
	for row in reader2:
		chr = row[0]
		istart = row[1]
		iend = row[2]
		hg19_intron = row[3]
		coverage = row[4]
		strand = row[5]
		GM12878_intron = chr + ":" + istart + strand + iend
		hg19_to_paternal[hg19_intron] = GM12878_intron
	
	for row in reader3:
		chr = row[0]
		istart = row[1]
		iend = row[2]
		hg19_intron = row[3]
		coverage = row[4]
		strand = row[5]
		GM12878_intron = chr + ":" + istart + strand + iend
		hg19_to_maternal[hg19_intron] = GM12878_intron	

	
	for row in reader1:
		intron = row[0]
		chr = row[1]
		strand = row[2]
		istart = int(row[3])
		iend = int(row[4])
		ilength = row[5]
		dn = row[6]
		dn_type = row[7]
		dn_type_score = row[8]
		bodymap_coverage = int(row[9])
		gm12878_coverage = int(row[10])
		hg19_cDNA_coverage = int(row[11])
		hg19_EST_coverage  = int(row[12])
		mm9_cDNA_coverage = int(row[13])
		mm9_EST_coverage = int(row[14])
		genecode_coverage = int(row[15])
		bodymap_seq = row[16]
		gm12878_seq = row[17]
		DR = row[18]
		intron_retention_exon = row[19].split(",")
		skipped_exons_names = row[20].split(",")
		alt_introns = row[21].split(",")
		alt_no_skipper_introns = row[22].split(",")
		alt_skipper_introns = row[23].split(",")
		alt_exon_variant_introns = row[24].split(",")
		shift = row[25].split(",")
		non_canonical_shift = row[26].split(",")
		
#		print intron_retention_exon, skipped_exons_names, alt_introns, alt_no_skipper_introns, alt_exon_variant_introns, shift

		retained_intron_fold = 0
		
		if row[19]!='NO':
		
			retained_istart_start = istart - 8
			retained_istart_end = istart + 8
			retained_iend_start = iend - 8  
			retained_iend_end = iend + 8
			
			retained_istart_seq_coverage = set([])
			retained_iend_seq_coverage = set([])
				
			for ali in hg19_bam.fetch(chr, retained_istart_start, retained_istart_end):
				start = ali.pos
				cigar = ali.cigar
				read_strand = "+"
				if ali.is_reverse:
					read_strand = "-"
				
				is_segment = False
				tstart = start
				for block in cigar:
					if block[0] == 0:
						tend = tstart + block[1]
						if tstart <= retained_istart_start and retained_istart_end <= tend and strand == read_strand:
							 is_segment = True
						tstart += block[1]
					elif block[0] == 2 or block[0] == 3 :
						tstart += block[1]
				
				if is_segment:
					retained_istart_seq_coverage.add(ali.seq)
				
			for ali in hg19_bam.fetch(chr, retained_iend_start, retained_iend_end):
				start = ali.pos
				cigar = ali.cigar
				read_strand = "+"
				if ali.is_reverse:
					read_strand = "-"
				
				is_segment = False
				tstart = start
				for block in cigar:
					if block[0] == 0:
						tend = tstart + block[1]
						if tstart <= retained_iend_start and retained_iend_end <= tend and strand == read_strand:
							 is_segment = True
						tstart += block[1]
					elif block[0] == 2 or block[0] == 3 :
						tstart += block[1]
				
				if is_segment:
					retained_iend_seq_coverage.add(ali.seq)
			
			retained_istart_coverage_hg19 = len(retained_istart_seq_coverage)
			retained_iend_coverage_hg19 = len(retained_iend_seq_coverage)
					
			retained_intron_coverage_hg19 = (Decimal(retained_istart_coverage_hg19 + retained_iend_coverage_hg19))/2
			

##### Se repite lo mismo para el paterno y materno de las GM12879, ademas se tubo que tomar en consideracion que la mayoria de estos reads son PE #####			

			retained_istart_seq_coverage_paternal = set([])
			retained_iend_seq_coverage_paternal = set([])

			try:

				paternal_intron = hg19_to_paternal[intron]
				paternal_chr = paternal_intron.split(":")[0]
						
#				paternal_istart = int(paternal_intron.split(":")[1].split(strand)[0])
				paternal_istart = int( re.split(":|\+|\-", paternal_intron)[1] )
				
				
#				paternal_iend = int(paternal_intron.split(":")[1].split(strand)[1])
				paternal_iend = int( re.split(":|\+|\-", paternal_intron)[2] )				
						
				paternal_retained_istart_start = paternal_istart - 8
				paternal_retained_istart_end = paternal_istart + 8
				paternal_retained_iend_start = paternal_iend - 8  
				paternal_retained_iend_end = paternal_iend + 8
						


							
				for ali in paternal_bam.fetch(paternal_chr, paternal_retained_istart_start, paternal_retained_istart_end):
					start = ali.pos
					cigar = ali.cigar
					read_strand = ""
					alignmet_strand = 1
					if ali.is_reverse:
						alignmet_strand = -1
					pair_orientation = 1
					if "/1" in ali.qname:   #En los datos los /2 son los forward y los /1 los reverse.
						pair_orientation = -1
							
					if alignmet_strand * pair_orientation == 1:
						read_strand = "+"
					elif alignmet_strand * pair_orientation == -1:
						read_strand = "-"
							
					is_segment = False
					tstart = start
					for block in cigar:
						if block[0] == 0:
							tend = tstart + block[1]
							if tstart <= paternal_retained_istart_start and paternal_retained_istart_end <= tend and strand == read_strand:
								 is_segment = True
							tstart += block[1]
						elif block[0] == 2 or block[0] == 3 :
							tstart += block[1]
							
					if is_segment:
						retained_istart_seq_coverage_paternal.add(ali.seq)


						
						
				for ali in paternal_bam.fetch(paternal_chr, paternal_retained_istart_start, paternal_retained_iend_end):
					start = ali.pos
					cigar = ali.cigar
					read_strand = ""
					alignmet_strand = 1
					if ali.is_reverse:
						alignmet_strand = -1
					pair_orientation = 1
					if "/1" in ali.qname:   #En los datos los /2 son los forward y los /1 los reverse.
						pair_orientation = -1
							
					if alignmet_strand * pair_orientation == 1:
						read_strand = "+"
					elif alignmet_strand * pair_orientation == -1:
						read_strand = "-"
							
					is_segment = False
					tstart = start
					for block in cigar:
						if block[0] == 0:
							tend = tstart + block[1]
							if tstart <= paternal_retained_istart_start and paternal_retained_iend_end <= tend and strand == read_strand:
								 is_segment = True
							tstart += block[1]
						elif block[0] == 2 or block[0] == 3 :
							tstart += block[1]
							
					if is_segment:
						retained_iend_seq_coverage_paternal.add(ali.seq)			


			except KeyError:
				pass

			retained_istart_seq_coverage_maternal = set([])
			retained_iend_seq_coverage_maternal = set([])

			try:

				maternal_intron = hg19_to_maternal[intron]
				maternal_chr = maternal_intron.split(":")[0]
						
#				maternal_istart = int(maternal_intron.split(":")[1].split(strand)[0])
				maternal_istart = int( re.split(":|\+|\-", maternal_intron)[1] )				
				
				
#				maternal_iend = int(maternal_intron.split(":")[1].split(strand)[1])
				maternal_iend = int( re.split(":|\+|\-", maternal_intron)[2] )

						
				maternal_retained_istart_start = maternal_istart - 8
				maternal_retained_istart_end = maternal_istart + 8
				maternal_retained_iend_start = maternal_iend - 8  
				maternal_retained_iend_end = maternal_iend + 8
						


							
				for ali in maternal_bam.fetch(maternal_chr, maternal_retained_istart_start, maternal_retained_istart_end):
					start = ali.pos
					cigar = ali.cigar
					read_strand = ""
					alignmet_strand = 1
					if ali.is_reverse:
						alignmet_strand = -1
					pair_orientation = 1
					if "/1" in ali.qname:   #En los datos los /2 son los forward y los /1 los reverse.
						pair_orientation = -1
							
					if alignmet_strand * pair_orientation == 1:
						read_strand = "+"
					elif alignmet_strand * pair_orientation == -1:
						read_strand = "-"
							
					is_segment = False
					tstart = start
					for block in cigar:
						if block[0] == 0:
							tend = tstart + block[1]
							if tstart <= maternal_retained_istart_start and maternal_retained_istart_end <= tend and strand == read_strand:
								 is_segment = True
							tstart += block[1]
						elif block[0] == 2 or block[0] == 3 :
							tstart += block[1]
							
					if is_segment:
						retained_istart_seq_coverage_maternal.add(ali.seq)
							
							
				for ali in maternal_bam.fetch(maternal_chr, maternal_retained_istart_start, maternal_retained_iend_end):
					start = ali.pos
					cigar = ali.cigar
					read_strand = ""
					alignmet_strand = 1
					if ali.is_reverse:
						alignmet_strand = -1
					pair_orientation = 1
					if "/1" in ali.qname:   #En los datos los /2 son los forward y los /1 los reverse.
						pair_orientation = -1
							
					if alignmet_strand * pair_orientation == 1:
						read_strand = "+"
					elif alignmet_strand * pair_orientation == -1:
						read_strand = "-"
							
					is_segment = False
					tstart = start
					for block in cigar:
						if block[0] == 0:
							tend = tstart + block[1]
							if tstart <= maternal_retained_istart_start and maternal_retained_iend_end <= tend and strand == read_strand:
								 is_segment = True
							tstart += block[1]
						elif block[0] == 2 or block[0] == 3 :
							tstart += block[1]
							
					if is_segment:
						retained_iend_seq_coverage_maternal.add(ali.seq)

			except KeyError:
				pass
				
				
			retained_istart_coverage_GM12878 = len(retained_istart_seq_coverage_paternal & retained_istart_seq_coverage_maternal)
			retained_iend_coverage_GM12878 = len(retained_iend_seq_coverage_paternal & retained_iend_seq_coverage_maternal)
						
			retained_intron_coverage_GM12878 = (Decimal(retained_istart_coverage_GM12878 + retained_iend_coverage_GM12878))/2
	

				
				
			retained_intron_fold = round(Decimal(retained_intron_coverage_hg19 + retained_intron_coverage_GM12878) / Decimal(bodymap_coverage + gm12878_coverage), 3)
			
		intron_retention_exon_coverage = []
		
		for i in intron_retention_exon:
			if retained_intron_fold != 0:
				intron_retention_exon_coverage.append(i + "|" + str(retained_intron_fold))
			

		if row[19] =='NO':	
			print "\t".join(row)
		
		else:
			if ",".join(intron_retention_exon_coverage)=="":	
				print "\t".join(row[:19]) + "\t" + "NO" + "\t" + "\t".join(row[20:])

			else:
				print "\t".join(row[:19]) + "\t" + ",".join(intron_retention_exon_coverage) + "\t" + "\t".join(row[20:])





				
				
				
			
#				print intron, hg19_to_maternal[intron], dn, round(Decimal(retained_intron_coverage)/Decimal(bodymap_coverage + gm12878_coverage), 3), intron_retention_exon, 
				
				
				

#python ~/my_src/Analisys/Alternative_splicing/intron_retention_coverage.py FINAL_TABLE.alternative_splicing  FINAL_TABLE.alternative_splicing.bed.paternal FINAL_TABLE.alternative_splicing.bed.maternal ../hg19/ALL/TOTAL.alignments.only-uniq.sort.bam ../GM12878/NA12878_Joel_Rozowsky/STRANDED/TOTAL_paternal/ALL/TOTAL.alignments.only-uniq.sort.bam ../GM12878/NA12878_Joel_Rozowsky/STRANDED/TOTAL_maternal/ALL/TOTAL.alignments.only-uniq.sort.bam 



#El problema aqui es que en el archivo de "retained_intron.bed.paternal" solo estan los intrones canonicos. Este archivo se ocupa para convertir las coordenadas de los intrones encontrados en el paterno y materno a hg19


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])

