import sys
import csv
from collections import defaultdict

def main (dn_maf, vertebrates_info, non_canonical_final_table):
	""" Hace un tacking de los dinucleotidos utilizando un maf """
	
	non_canonical_introns = []
	
	for row in csv.reader(open(non_canonical_final_table), delimiter = '\t'):
		
		name = row[0]  	
		intron = row[1]
		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
		ilength = row[6]
		dn = row[7]
		dn_type = row[8]
		dn_type_score = float(row[9])
		bodymap_coverage = int(row[10])
		gm12878_coverage = int(row[11])
		hg19_cDNA_coverage = int(row[12])
		hg19_EST_coverage  = int(row[13])
		mm9_cDNA_coverage = int(row[14])
		mm9_EST_coverage = int(row[15])
		genecode_coverage = int(row[16])
		TOTAL_tissues_coverage = int(row[17])
		n_tissues = row[18]	
		intron_retention_exon = row[19].split(",")
		skipped_exons_names = row[20].split(",")
		alt_introns = row[21].split(",")
		alt_no_skipper_introns = row[22].split(",")
		alt_skipper_introns = row[23].split(",")
		alt_exon_variant_introns = row[24].split(",")
		shift = row[25].split(",")
		non_canonical_shift = row[26].split(",")
		
		non_canonical_introns.append(" ".join(row[:10]))
	
	genomes = []
	
	for row in csv.reader(open(vertebrates_info), delimiter = ' '):
		genome = row[-1]
		specie = (" ").join(row[:-3])
		
		genomes.append(genome)
	
	dn_pos = 0
	maf_block = -1
	
	dn_genomes = defaultdict(list)	
	
		
	for row in csv.reader(open(dn_maf), delimiter = ' '):
		
		species_key = set([])
		
		if len(row)!=0 and row[0]=="s":
			
			while True:
			  try:
				row.remove("")
			  except ValueError:
				break
			
			genome = row[1].split(".")[0]
			nt = row[-1].upper()
			
			if genome == "hg19":
				dn_pos +=1

			if dn_pos > 4:
				dn_pos = 1
				maf_block += 1
				
				print non_canonical_introns[maf_block]
				
				for g in genomes:
					
					
					nt_dn_pos_dict = {1:"-", 2:"-", 3:"-", 4:"-"}
					
					for n in dn_genomes[g]:
						nt_dn_pos = int(n[-1])
						nt = n[:-1]
						nt_dn_pos_dict[nt_dn_pos] = nt					
					
					print g, nt_dn_pos_dict[1], nt_dn_pos_dict[2], nt_dn_pos_dict[3], nt_dn_pos_dict[4] 
				
				print ""	
				dn_genomes = defaultdict(list) 
				
			dn_genomes[genome].append(nt + str(dn_pos))

	print non_canonical_introns[-1]
	
	for g in genomes:
					
		nt_dn_pos_dict = {1:"-", 2:"-", 3:"-", 4:"-"}
					
		for n in dn_genomes[g]:
			nt_dn_pos = int(n[-1])
			nt = n[:-1]
			nt_dn_pos_dict[nt_dn_pos] = nt					
					
		print g, nt_dn_pos_dict[1], nt_dn_pos_dict[2], nt_dn_pos_dict[3], nt_dn_pos_dict[4] 

	#print dn_genomes
							

		
		
			
			
		
		

	
if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])	
