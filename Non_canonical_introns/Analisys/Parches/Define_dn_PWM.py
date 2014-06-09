import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable = []
		
def Genomictabulator(fasta):
	
	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

	f.close() 

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0		

def PWM_to_dict(file):
	reader = csv.reader(open(file), delimiter = '\t')
	header = reader.next()
	header_dict = {}
	col = 0
	
	matrix = {}
	
	for name in header:
		header_dict[name] = col
		col += 1
	
	A_frec = []
	C_frec = []
	G_frec = []
	T_frec = []
	
	for row in reader:
		A = row[header_dict["A"]]
		C = row[header_dict["C"]]
		G = row[header_dict["G"]]
		T = row[header_dict["T"]]
		
		A_frec.append(float(A))
		C_frec.append(float(C))
		G_frec.append(float(G))
		T_frec.append(float(T))
			
	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	
	return matrix

def DR_SJ_slicer(intron, ichr, strand, istart, iend, Genome):
	
	'''Esta funcion busca intrones canonicos (GTAG, ATAC, GCAG) al correr el intron por sus directos repetidos''' 
	
	U2_5_DR = []
	U2_3_DR = []
	U12_5_DR = []
	U12_3_DR = []
	istart_DR = [istart]
	iend_DR = [iend]
			
	L = 100         #Solo permite que se corra L pares de bases para buscar DR 		

	#Extrayendo regiones exonicas colindantes
				
	SJ5U = Genome[ichr][istart-L : istart].lower()
	SJ5D = Genome[ichr][istart : istart+L].lower()
	SJ3U = Genome[ichr][iend-L : iend].lower()
	SJ3D = Genome[ichr][iend : iend+L].lower()
	
	U2_5 = Genome[ichr][(istart-2):(istart+5)]
	U2_3 = Genome[ichr][(iend-7):(iend+1)]
	U12_5 = Genome[ichr][(istart-1):(istart+8)]
	U12_3 = Genome[ichr][(iend-4):(iend+2)]
	
	if strand == "-":
		SJ5U = Genome[ichr][iend : iend+L].lower().reverse_complement() 
		SJ5D = Genome[ichr][iend-L : iend].lower().reverse_complement()
		SJ3U = Genome[ichr][istart : istart+L].lower().reverse_complement()
		SJ3D = Genome[ichr][istart-L : istart].lower().reverse_complement()
		
		
		U2_5 = Genome[ichr][(iend-5):(iend+2)]
		U2_3 = Genome[ichr][(istart-1):(istart+7)]
		U12_5 = Genome[ichr][(iend-8):(iend+1)]
		U12_3 = Genome[ichr][(istart-2):(istart+4)]
		
		U2_5 = U2_5.reverse_complement()
		U2_3 = U2_3.reverse_complement()
		U12_5 = U12_5.reverse_complement()
		U12_3 = U12_3.reverse_complement()
		
	U2_5 = str(U2_5).upper()
	U2_3 = str(U2_3).upper()
	U12_5 = str(U12_5).upper()
	U12_3 = str(U12_3).upper()
	
	U2_5_DR.append(U2_5)
	U2_3_DR.append(U2_3)
	U12_5_DR.append(U12_5)
	U12_3_DR.append(U12_3)			
	
	DRU = 0
	DRD = 0
								
	#Contando directos repetidos y generando intrones no consenso alternativos
		
	try:
		while SJ5U[L-1-DRU]==SJ3U[L-1-DRU]:
			DRU += 1
			if strand == "+":
				U2_5 = Genome[ichr][(istart-DRU-2):(istart+5-DRU)]
				U2_3 = Genome[ichr][(iend-7-DRU):(iend-DRU+1)]
				U12_5 = Genome[ichr][(istart-DRU-1):(istart+8-DRU)]
				U12_3 = Genome[ichr][(iend-4-DRU):(iend-DRU+2)]
								
				U2_5 = str(U2_5).upper()
				U2_3 = str(U2_3).upper()
				U12_5 = str(U12_5).upper()
				U12_3 = str(U12_3).upper()							
				
				U2_5_DR.append(U2_5)
				U2_3_DR.append(U2_3)
				U12_5_DR.append(U12_5)
				U12_3_DR.append(U12_3)
				istart_DR.append(istart-DRU)
				iend_DR.append(iend-DRU)
				
			elif strand == "-":
				U2_5 = Genome[ichr][(iend-5+DRU):(iend+2+DRU)]
				U2_3 = Genome[ichr][(istart-1+DRU):(istart+7+DRU)]
				U12_5 = Genome[ichr][(iend-8+DRU):(iend+1+DRU)]
				U12_3 = Genome[ichr][(istart-2+DRU):(istart+4+DRU)]
				
				U2_5 = str(U2_5.reverse_complement()).upper()
				U2_3 = str(U2_3.reverse_complement()).upper()
				U12_5 = str(U12_5.reverse_complement()).upper()
				U12_3 = str(U12_3.reverse_complement()).upper()					
				
				U2_5_DR.append(U2_5)
				U2_3_DR.append(U2_3)
				U12_5_DR.append(U12_5)
				U12_3_DR.append(U12_3)
				istart_DR.append(istart+DRU)
				iend_DR.append(iend+DRU)								

			if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
				break
	except IndexError:
		pass 
	try:
		while SJ5D[DRD]==SJ3D[DRD]:
			DRD += 1
			if strand == "+":
				U2_5 = Genome[ichr][(istart+DRD-2):(istart+5+DRD)]
				U2_3 = Genome[ichr][(iend-7+DRD):(iend+DRD+1)]
				U12_5 = Genome[ichr][(istart+DRD-1):(istart+8+DRD)]
				U12_3 = Genome[ichr][(iend-4+DRD):(iend+DRD+2)]
								
				U2_5 = str(U2_5).upper()
				U2_3 = str(U2_3).upper()
				U12_5 = str(U12_5).upper()
				U12_3 = str(U12_3).upper()							
				
				U2_5_DR.append(U2_5)
				U2_3_DR.append(U2_3)
				U12_5_DR.append(U12_5)
				U12_3_DR.append(U12_3)
				istart_DR.append(istart+DRD)
				iend_DR.append(iend+DRD)
				
			elif strand == "-":
				U2_5 = Genome[ichr][(iend-5-DRD):(iend+2-DRD)]
				U2_3 = Genome[ichr][(istart-1-DRD):(istart+7-DRD)]
				U12_5 = Genome[ichr][(iend-8-DRD):(iend+1-DRD)]
				U12_3 = Genome[ichr][(istart-2-DRD):(istart+4-DRD)]
				
				U2_5 = str(U2_5.reverse_complement()).upper()
				U2_3 = str(U2_3.reverse_complement()).upper()
				U12_5 = str(U12_5.reverse_complement()).upper()
				U12_3 = str(U12_3.reverse_complement()).upper()					
				
				U2_5_DR.append(U2_5)
				U2_3_DR.append(U2_3)
				U12_5_DR.append(U12_5)
				U12_3_DR.append(U12_3)
				istart_DR.append(istart-DRD)
				iend_DR.append(iend-DRD)
								
			if SJ5D[DRD]!=SJ3D[DRD]:
				break
	except IndexError:
		pass
	
	return (U2_5_DR, U2_3_DR, U12_5_DR, U12_3_DR, istart_DR, iend_DR)
		
	
		
def main(introns_final_table, U2_GTAG_5_file, U2_GTAG_3_file, U2_GCAG_5_file, U2_GCAG_3_file, U12_GTAG_5_file, U12_GTAG_3_file, U12_ATAC_5_file, U12_ATAC_3_file  ):
	Genome = dict(SeqTable)

	csv.field_size_limit(1000000000)

	reader = csv.reader(open(introns_final_table), delimiter = ' ')
	

# Ventanas a utilizar	
#  U2 = EE IIIII ------ IIIIIII E
#       2 	 5             7    1
# U12 = E IIIIIIII ------ IIII EE
#       1     8             4   2

	
	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)
	U2_GCAG_5 = PWM_to_dict(U2_GCAG_5_file)
	U2_GCAG_3 = PWM_to_dict(U2_GCAG_3_file)
	U12_GTAG_5 = PWM_to_dict(U12_GTAG_5_file)
	U12_GTAG_3 = PWM_to_dict(U12_GTAG_3_file)
	U12_ATAC_5 = PWM_to_dict(U12_ATAC_5_file)
	U12_ATAC_3 = PWM_to_dict(U12_ATAC_3_file)

		
	U2_GTAG_5_win = range(7)
	U2_GTAG_5_max_score = 0
	index = 1
	for N in U2_GTAG_5_win:
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])
		index += 1
	
	U2_GTAG_3_win = range(8)
	U2_GTAG_3_max_score = 0
	index = 7
	for N in U2_GTAG_3_win:
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])
		index += 1
	
	U2_GCAG_5_win = range(7)
	U2_GCAG_5_max_score = 0
	index = 1
	for N in U2_GCAG_5_win:
		U2_GCAG_5_max_score += max(U2_GCAG_5['A'][index], U2_GCAG_5['C'][index], U2_GCAG_5['T'][index], U2_GCAG_5['G'][index])
		index += 1
	
	U2_GCAG_3_win = range(8)
	U2_GCAG_3_max_score = 0
	index = 7
	for N in U2_GCAG_3_win:
		U2_GCAG_3_max_score += max(U2_GCAG_3['A'][index], U2_GCAG_3['C'][index], U2_GCAG_3['T'][index], U2_GCAG_3['G'][index])
		index += 1

	U12_GTAG_5_win = range(9)
	U12_GTAG_5_max_score = 0
	index = 2
	for N in U12_GTAG_5_win:
		U12_GTAG_5_max_score += max(U12_GTAG_5['A'][index], U12_GTAG_5['C'][index], U12_GTAG_5['T'][index], U12_GTAG_5['G'][index])
		index += 1
	
	U12_GTAG_3_win = range(6)
	U12_GTAG_3_max_score = 0
	index = 10
	for N in U12_GTAG_3_win:
		U12_GTAG_3_max_score += max(U12_GTAG_3['A'][index], U12_GTAG_3['C'][index], U12_GTAG_3['T'][index], U12_GTAG_3['G'][index])
		index += 1
	
	U12_ATAC_5_win = range(9)
	U12_ATAC_5_max_score = 0
	index = 2
	for N in U12_ATAC_5_win:
		U12_ATAC_5_max_score += max(U12_ATAC_5['A'][index], U12_ATAC_5['C'][index], U12_ATAC_5['T'][index], U12_ATAC_5['G'][index])
		index += 1
	
	U12_ATAC_3_win = range(6)
	U12_ATAC_3_max_score = 0
	index = 10

	for N in U12_ATAC_3_win:
		U12_ATAC_3_max_score += max(U12_ATAC_3['A'][index], U12_ATAC_3['C'][index], U12_ATAC_3['T'][index], U12_ATAC_3['G'][index])
		index += 1	
		
	
#	print U2_GTAG_5_max_score, U2_GTAG_3_max_score, U2_GCAG_5_max_score, U2_GCAG_3_max_score, U12_GTAG_5_max_score, U12_GTAG_3_max_score, U12_ATAC_5_max_score, U12_ATAC_3_max_score


	for row in reader:
		
#		intron = row[0]

		chr = row[2]
		strand = row[3]
		istart = int(row[4])
		iend = int(row[5])
#		ilength = row[5]
		dn = row[7]
		dn_type = row[8]
		dn_type_score = row[9]
#		bodymap_coverage = int(row[9])
#		gm12878_coverage = int(row[10])
#		hg19_cDNA_coverage = int(row[11])
#		hg19_EST_coverage  = int(row[12])
#		mm9_cDNA_coverage = int(row[13])
#		mm9_EST_coverage = int(row[14])
#		genecode_coverage = int(row[15])
#		bodymap_seq = row[16]
#		gm12878_seq = row[17]

		

			 
		istart = int(istart)
		iend = int(iend)

		u2_5 = Genome[chr][(istart-2):(istart+5)]
		u2_3 = Genome[chr][(iend-7):(iend+1)]
		u12_5 = Genome[chr][(istart-1):(istart+8)]
		u12_3 = Genome[chr][(iend-4):(iend+2)]
								
		u2_5 = str(u2_5).upper()
		u2_3 = str(u2_3).upper()
		u12_5 = str(u12_5).upper()
		u12_3 = str(u12_3).upper()							

				
		if strand == "-":
			u2_5 = Genome[chr][(iend-5):(iend+2)]
			u2_3 = Genome[chr][(istart-1):(istart+7)]
			u12_5 = Genome[chr][(iend-8):(iend+1)]
			u12_3 = Genome[chr][(istart-2):(istart+4)]
				
			u2_5 = str(u2_5.reverse_complement()).upper()
			u2_3 = str(u2_3.reverse_complement()).upper()
			u12_5 = str(u12_5.reverse_complement()).upper()
			u12_3 = str(u12_3.reverse_complement()).upper()		
			

		U2_GTAG_5_score = 0
		index = 1 
		for N in u2_5:
			try:
				frec = U2_GTAG_5[N][index]
				U2_GTAG_5_score += frec
				index += 1
			except KeyError:
				pass

		U2_GTAG_3_score = 0
		index = 7 
		for N in u2_3:
			try:
				frec = U2_GTAG_3[N][index]
				U2_GTAG_3_score += frec
				index += 1					 
			except KeyError:
				pass
				
		U12_GTAG_5_score = 0
		index = 2 
		for N in u12_5:
			try:
				frec = U12_GTAG_5[N][index]
				U12_GTAG_5_score += frec
				index += 1
			except KeyError:
				pass	

		U12_GTAG_3_score = 0
		index = 10 
		for N in u12_3:
			try:
				frec = U12_GTAG_3[N][index]
				U12_GTAG_3_score += frec
				index += 1						
			except KeyError:
				pass
				
		U2_GTAG_total_score = percent( (U2_GTAG_5_score + U2_GTAG_3_score), (U2_GTAG_5_max_score + U2_GTAG_3_max_score) )
		U12_GTAG_total_score = percent( (U12_GTAG_5_score + U12_GTAG_3_score), (U12_GTAG_5_max_score + U12_GTAG_3_max_score) )
			
		dn_type = "U2_GTAG"
		score = U2_GTAG_total_score
		if U12_GTAG_total_score > U2_GTAG_total_score:
			dn_type = "U12_GTAG"
			score = U12_GTAG_total_score
		
		if row[8]=="U12_ATAC":
			print " ".join(row)
		else:
		

#			print " ".join(row[:7]), dn_type, score, " ".join(row[9:])
			print " ".join(row[:8]), dn_type, score, " ".join(row[10:])





#python ~/my_src/Analisys/Parches/Define_dn_PWM.py ~/db/genome/hg19.fa introns.final_table.hg19.fixed.tags /home/geparada/db/PWM/hg19_GT_AG_U2_5.good.matrix /home/geparada/db/PWM/hg19_GT_AG_U2_3.good.matrix /home/geparada/db/PWM/hg19_GC_AG_U2_5.good.matrix /home/geparada/db/PWM/hg19_GC_AG_U2_3.good.matrix /home/geparada/db/PWM/hg19_GT_AG_U12_5.good.matrix /home/geparada/db/PWM/hg19_GT_AG_U12_3.good.matrix /home/geparada/db/PWM/hg19_AT_AC_U12_5.good.matrix /home/geparada/db/PWM/hg19_AT_AC_U12_3.good.matrix		

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
