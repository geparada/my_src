import sys
import csv
from collections import defaultdict
from decimal import *

#Esta version tiene una modificacion para que la variantes de splicing que ocurran en el ultiumo exon sean bien clasificadas.



def main(final_table, genecode_bed):
	
	getcontext().prec = 20
		
	istart_coverage = defaultdict(int)
	iend_coverage = defaultdict(int)	
	alt_start_end = defaultdict(list)
	
	intron_skipped_exons = {}
	intron_retention_exons = {}

	non_canonical = []
	
	TOTAL_intron = defaultdict(list)
	TOTAL_intron_set = set([])
	
	exons = set([])
	exons_coverage = defaultdict(list)
	
	introns_genecode = []	
	
	
	for row in csv.reader(open(genecode_bed), delimiter = '\t'):

		qstarts = map (int, row[11].strip(",").split(","))                      
		blocksizes = map(int, row[10].strip(",").split(","))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chr = row[0]
		
		blocks_number = len(blocksizes)
		blockcount = 0

		for q, b in zip(qstarts, blocksizes):
			exon_pos = "middle"
			blockcount += 1
			if blockcount == 1:
				exon_pos = "first"
			if blockcount == blocks_number:
				exon_pos = "last"
			
			estart = start+q
			eend = start+q+b
			exon = chr + ":" + str(estart) + strand + str(eend)

			exons.add((exon, chr, strand, estart, eend, exon_pos))
	

	for row in csv.reader(open(final_table), delimiter = ' '):

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
		tissues_coverage = int(row[16])
		n_tissues = row[17]
		tissues = row[18]
	

		istart_coverage[chr+strand+str(istart)] += bodymap_coverage + gm12878_coverage + tissues_coverage
		iend_coverage[chr+strand+str(iend)] += bodymap_coverage	+ gm12878_coverage  + tissues_coverage
		
		alt_start_end[chr+strand+str(istart)+"s"].append(row)
		alt_start_end[chr+strand+str(iend)+"e"].append(row)
		
		
#		if dn !="GTAG" and dn != "GCAG" and dn != "ATAC":
		if 0==0:
			
			non_canonical.append(row)
#		else:
		TOTAL_intron[chr+strand].append(row)
		TOTAL_intron_set.add(intron)
			

	
	for e in exons:
		exon = e[0]
		chr = e[1]
		strand = e[2]
		estart = e[3]
		eend = e[4]
		exon_pos = e[5]	
		coverages = []
		coverages.append(iend_coverage[chr+strand+str(estart)])
		coverages.append(istart_coverage[chr+strand+str(eend)])

		
		if 0 in coverages:
			if exon_pos=="first" or exon_pos=="last": 
				coverages.remove(0)
			else:
				coverages = [0]
		exon_coverage = float(sum(coverages))/float(len(coverages))
		exons_coverage[chr+strand].append((exon, chr, strand, estart, eend, exon_pos, exon_coverage))
		
		
	for row in non_canonical:
		
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
		tissues_coverage = int(row[16])
		n_tissues = row[17]
		tissues = row[18]
		
		skipped_exons = []
		retention_exons = []
		for e in exons_coverage[chr+strand]:
			exon = e[0]
			chr = e[1]
			strand = e[2]
			estart = e[3]
			eend = e[4]
			exon_pos = e[5]	
			exon_coverage = e[6]
#			if exon_pos == "middle":   #Aqui me equiboque!!, los last tambien es mejor considerarlos como exon skiped!
			if istart<=estart and eend<=iend:
				skipped_exons.append((chr, strand, estart, eend, str(round(Decimal(exon_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) ))
			if estart<=istart and iend<=eend:
				retention_exons.append((chr, strand, estart, eend))
		
		intron_retention_exons[intron] = retention_exons								
		intron_skipped_exons[intron] = skipped_exons
		

	
	for row in non_canonical:
		
		
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
		tissues_coverage = int(row[16])
		n_tissues = row[17]
		tissues = row[18]
		
		info_start = []
		
		if alt_start_end.has_key(chr+strand+str(istart)+"s"):
			info_start = alt_start_end[chr+strand+str(istart)+"s"]
		
		tuple_info_start = map(lambda x: tuple(x), info_start)
		
		set_info_start = set(tuple_info_start)
		
		info_end = []
		
				
		if alt_start_end.has_key(chr+strand+str(iend)+"e"):
			info_end = alt_start_end[chr+strand+str(iend)+"e"]
				
		tuple_info_end = map(lambda x: tuple(x), info_end)
		
		set_info_end = set(tuple_info_end)
		
		alt_info = (set_info_start | set_info_end) - set([tuple(row)])

		
		alt_introns = []
		alt_skipper_introns = []
		alt_no_skipper_introns = []
		alt_exon_variant_introns = []
		
				
		alt_canonical = False    #indicador para saber si hay un intron canonico con el donor o aceptor del no-canonico
					
		for alt in alt_info:
			alt_intron = alt[0]
			alt_chr = alt[1]
			alt_strand = alt[2]
			alt_istart = int(alt[3])
			alt_iend = int(alt[4])
			alt_dn = alt[6]
			alt_bodymap_coverage = int(alt[9])
			alt_gm12878_coverage = int(alt[10])
			alt_tissues_coverage = int(alt[16])			
			skipped_exons = intron_skipped_exons[intron]
			
			if alt_dn == "GTAG" or alt_dn == "GCAG" or alt_dn == "ATAC" :
				alt_canonical = True
			
			alt_skipped_exons = []
				
			for e in exons_coverage[chr+strand]:
				exon = e[0]
				chr = e[1]
				strand = e[2]
				estart = e[3]
				eend = e[4]
				exon_pos = e[5]	
				exon_coverage = e[6]
#				if exon_pos == "middle":
				if alt_istart<=estart and alt_iend>=eend:
					alt_skipped_exons.append((chr, strand, estart, eend, str(round(Decimal(exon_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) ))

			is_alt_no_skipper_intron = False			
			for e in skipped_exons:           #Modulo para encontrar las variantes de intrones que son resultado no saltarse un intron
				if (alt_istart in e) or (alt_iend in e):
					is_alt_no_skipper_intron = True
			if is_alt_no_skipper_intron == True and (alt_istart != istart or alt_iend != iend) :
				if strand == "+":
					alt_no_skipper_introns.append(alt_intron + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(alt_istart - istart) + "|" + str(alt_iend - iend))
				if strand == "-":
					alt_no_skipper_introns.append(str(alt_intron) + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(iend - alt_iend) + "|" + str(istart - alt_istart))	


			if set(alt_skipped_exons) == set(skipped_exons) and (alt_istart != istart or alt_iend != iend):         #Modulo para en contrar variantes de intrones que son producidas por saltarse mas exones que el no canonico
				if strand == "+":
					alt_exon_variant_introns.append(alt_intron + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(alt_istart - istart) + "|" + str(alt_iend - iend))
				if strand == "-":
					alt_exon_variant_introns.append(str(alt_intron) + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage  + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(iend - alt_iend) + "|" + str(istart - alt_istart))	

															
			if len(set(alt_skipped_exons) -set(skipped_exons))!=0 and (alt_istart != istart or alt_iend != iend):         #Modulo para en contrar variantes de intrones que son producidas por saltarse mas exones que el no canonico
				if strand == "+":
					alt_skipper_introns.append(alt_intron + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(alt_istart - istart) + "|" + str(alt_iend - iend))
				if strand == "-":
					alt_skipper_introns.append(str(alt_intron) + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(iend - alt_iend) + "|" + str(istart - alt_istart))																	

							
			if strand == "+" and (alt_istart != istart or alt_iend != iend):           #Modulo para encontrar el total de las variates 5' y 3'
				alt_introns.append(alt_intron + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(alt_istart - istart) + "|" + str(alt_iend - iend))
			if strand == "-" and (alt_istart != istart or alt_iend != iend):
				alt_introns.append(str(alt_intron) + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage  + tissues_coverage), 3)) + "|" + str(iend - alt_iend) + "|" + str(istart - alt_istart))



		shift = []
		non_canonical_shift = []
		retention_exons = []
		retention_exons_names = []
		skipped_exons = []
		skipped_exons_names = []
		
		try:
			retention_exons = intron_retention_exons[intron]

		except KeyError:
			pass												
			
		for e in retention_exons:
			chr = e[0]
			strand = e[1]
			estart = e[2]
			eend = e[3]
			name = chr + ":" + str(estart) + strand + str(eend)
			retention_exons_names.append(name)
			
		try:
			skipped_exons = intron_skipped_exons[intron]

		except KeyError:
			pass												
			
		for e in skipped_exons:
			chr = e[0]
			strand = e[1]
			estart = e[2]
			eend = e[3]
			fold = e[4]
			name = chr + ":" + str(estart) + strand + str(eend)
			skipped_exons_names.append((name + "|" + fold)) 
					
		
		if alt_canonical == False:          #solo busca si hay shifts cuando no se encontro un intron canonico alternativo
			TOTAL_intron_info = TOTAL_intron[chr+strand]
			for alt in TOTAL_intron_info:

				alt_intron = alt[0]
				alt_chr = alt[1]
				alt_strand = alt[2]
				alt_istart = int(alt[3])
				alt_iend = int(alt[4])
				alt_dn = alt[6]
				alt_bodymap_coverage = int(alt[9])
				alt_gm12878_coverage = int(alt[10])
					
				if (alt_istart-10) <= istart and istart <= (alt_istart+10) and (alt_iend-10) <= iend and iend <= (alt_iend + 10):
					if (alt_dn == "GTAG" or alt_dn == "GCAG" or alt_dn == "ATAC") and (alt_istart != istart or alt_iend != iend) :					
						if strand== "+":
							shift.append(alt_intron + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage + tissues_coverage), 3)) + "|" + str(alt_istart - istart) + "|" + str(alt_iend - iend))
						if strand== "-":
							shift.append(str(alt_intron) + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage + tissues_coverage), 3)) + "|" + str(iend - alt_iend) + "|" + str(istart - alt_istart))

					elif (alt_istart - istart) !=0 and (alt_iend - iend) != 0 and (alt_istart != istart or alt_iend != iend):   #Para evitar que se reporte el intron no canonico como shift de si mismo, o que alt 3 o 5
						if strand== "+":
							non_canonical_shift.append(alt_intron + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage + tissues_coverage), 3)) + "|" + str(alt_istart - istart) + "|" + str(alt_iend - iend))
						if strand== "-":
							non_canonical_shift.append(str(alt_intron) + "|" + alt_dn + "|" + str(round(Decimal(alt_bodymap_coverage + alt_gm12878_coverage + alt_tissues_coverage)/Decimal(bodymap_coverage + gm12878_coverage + tissues_coverage), 3)) + "|" + str(iend - alt_iend) + "|" + str(istart - alt_istart))
						

		if bodymap_coverage + gm12878_coverage == 0:
			retention_exons_names = ""
		
		alt_table = [ ','.join(retention_exons_names), ','.join(skipped_exons_names), ','.join(alt_introns), ','.join(alt_no_skipper_introns), ','.join(alt_skipper_introns),  ','.join(alt_exon_variant_introns), ','.join(shift), ','.join(non_canonical_shift)]
				
		print "\t".join(row) + "\t" + "\t".join(["NO" if x=="" else x for x in alt_table])
						
						

	



#python ~/my_src/Analisys/Alternative_splicing/alternative_splicing.py FINAL_TABLE ~/db/transcriptome/hg19/Gene_models/gencode/v11/gencode.v11.annotation.bed12



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
