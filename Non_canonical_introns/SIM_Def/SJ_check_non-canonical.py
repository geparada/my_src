import sys
import csv
import random 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable=[]


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ..."	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable.append(table)

	f.close()  



def SJcheck(bed12):
	'''Verifica si es que los reads simulados de SJ alinearon correctamente
	*Input file = SJ.bed12
	*Usage = python SJ_check file
	*Output =  '''  

	print >> sys.stderr, "Revisando BED12"     

  	Genome = dict(SeqTable)

        

	reader1 = csv.reader(open(bed12), dialect='excel-tab' )
         
		# row[0]  chr
		# row[1]  alignment start
		# row[2]  alignment end
		# row[3]  name (Refseq=intron0<>intron1<>intron2=anchor0<>anchor1<>anchor2)
		# row[4]    
		# row[5]  strand
		# row[6]  aligment start
		# row[7]  aligment end
		# row[8]  
		# row[9] blocknum
		# row[10] blocksizes
		# row[11] qstarts         
    
    
   
	for row in reader1:

		list_anchors_finded = []
		list_anchors_sim = []
        
		q = map (int, row[11].strip(",").split(","))                      # obteniendo lista int de las distintas variables de la linea psl
		b = map(int, row[10].strip(",").split(","))
		start = int(row[1])
		bn = int(row[9])	                  #El blocknum es importante para saber cuantos intrones hay en el alineamiento
		chr = row[0]	
	                                                   
		intron1 = ''
		intron2 = ''
		intron3 = ''

		anchor1 = ''		
		anchor2 = ''	
		anchor3	= ''	
		

		if bn>=2 and ((start + q[1]) - (start + q[0] + b[0]))>=40:
			istart = start + q[0] + b[0]
			iend = start + q[1]
			strand = row[5]
			intron1 = row[0] + ':' +  str(istart) + strand + str(iend)
			anchor1 = str(min([b[0], b[1]]))
			list_anchors_finded.append((intron1,anchor1))

   
		if bn>=3 and ((start + q[2]) - (start + q[1] + b[1]))>=40:
			istart = start + q[1] + b[1]
			iend = start + q[2]
			strand = row[5]
			intron2 = row[0] + ':' +  str(istart) + strand + str(iend)
			anchor2 = str(min([b[1], b[2]]))
			list_anchors_finded.append((intron2,anchor2))


		if bn>=4 and ((start + q[3]) - (start + q[2] + b[2]))>=40:
			istart = start + q[2] + b[2]
			iend = start + q[3]
			strand = row[5]
			intron3 = row[0] + ':' +  str(istart) + strand + str(iend)
			anchor3 = str(min([b[2], b[3]]))
			list_anchors_finded.append((intron3,anchor3))

        
		introns_finded = [intron1, intron2, intron3]
		introns_finded_DR = [intron1, intron2, intron3]        

		introns_sim = row[3].split('=')[1].split('<>')
		anchors_sim = row[3].split('=')[2].split('<>')

		for i, a in zip(introns_sim, anchors_sim):
			if i!='':
				min_a = min(map(int,a.split("-")))
				list_anchors_sim.append((i,min_a))                                   
	
		bn_sim = len(introns_sim) + 1

		read_length = sum(map(int, anchors_sim[0].split('-')))
	
		if bn_sim==3:
			read_length = read_length + int(anchors_sim[1].split('-')[1])
		if bn_sim==4:
			read_length = read_length + int(anchors_sim[1].split('-')[1]) + int(anchors_sim[2].split('-')[1]) 	

#		print row
		DR_list = []
		DRU_list = []
		DRD_list = []		


		for i in introns_finded:
			ichr = i.split(':')[0]		
			istart = 0
			iend = 0
			strand = ''
		
			if '+' in i:
				istart = int(i.split(':')[1].split('+')[0])
				iend = int(i.split(':')[1].split('+')[1])
				strand = '+'
			elif '-' in i:
				istart = int(i.split(':')[1].split('-')[0])
				iend = int(i.split(':')[1].split('-')[1])
				strand = '-'

			L = 100         #Solo permite que se corra L pares de bases para buscar DR 		

			if ichr != '':

				SJ5U = Genome[ichr][istart-L : istart].lower()
				SJ5D = Genome[ichr][istart : istart+L].lower()
				SJ3U = Genome[ichr][iend-L : iend].lower()
				SJ3D = Genome[ichr][iend : iend+L].lower()
			
				if '-' in i:
					SJ5U = Genome[ichr][iend : iend+L].lower().reverse_complement() 
					SJ5D = Genome[ichr][iend-L : iend].lower().reverse_complement()
					SJ3U = Genome[ichr][istart : istart+L].lower().reverse_complement()
					SJ3D = Genome[ichr][istart-L : istart].lower().reverse_complement()
			
				DRU = 0
				DRD = 0
										
				
				try:
					while SJ5U[L-1-DRU]==SJ3U[L-1-DRU]:
						DRU = (DRU + 1)
						if strand == "+":
							introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart-DRU) + strand + str(iend-DRU)]
						elif strand == "-":
							introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart+DRU) + strand + str(iend+DRU)]
						if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
							break
				except IndexError:
					pass 

				try:
					while SJ5D[DRD]==SJ3D[DRD]:
						DRD = (DRD + 1)
						if strand == "+":
							introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart+DRD) + strand + str(iend+DRD)]
						elif strand == "-":
							introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart-DRD) + strand + str(iend-DRD)]
						if SJ5D[DRD]!=SJ3D[DRD]:
							break
				except IndexError:
					pass

				DR_list.append(DRU+DRD)
				DRU_list.append(DRU)
				DRD_list.append(DRD)



		#Definiendo conjuntos de respuestas generadas apartir de los directos repetidos para cada intron encontrado

		intron1_DR = [introns_finded_DR[0]]        
		intron2_DR = [introns_finded_DR[1]]
		intron3_DR = [introns_finded_DR[2]]

		try:
			if DR_list[0] != 0:
				intron1_DR = intron1_DR + introns_finded_DR[3:3+DR_list[0]]
		except IndexError:
			pass 		

		try:
			if DR_list[1] != 0:
				intron2_DR = intron2_DR + introns_finded_DR[3+DR_list[0]:3+DR_list[0]+DR_list[1]]
		except IndexError:
			pass 	

		try:
			if DR_list[2] != 0:
				intron3_DR = intron3_DR + introns_finded_DR[3+DR_list[0]+DR_list[1]:3+DR_list[0]+DR_list[1]+DR_list[2]]
		except IndexError:
			pass
				
		 			
		OK = set(introns_sim) & set(introns_finded_DR)
		Wrong = set(introns_finded_DR[:3]) - set([""])
		OK_finded = []
		

		for i in OK:
			if i in intron1_DR:
				Wrong = Wrong - set(intron1_DR)
				OK_finded += [introns_finded_DR[0]]

			if i in intron2_DR:
				Wrong = Wrong - set(intron2_DR)
				OK_finded += [introns_finded_DR[1]]

			if i in intron3_DR:
				Wrong = Wrong - set(intron3_DR)
				OK_finded += [introns_finded_DR[2]] 

#		if len(introns_sim)==3:
#			print introns_sim, introns_finded_DR[:3], "OK", OK, "Wrong", Wrong   #, DR_list, intron2_DR 

		N_intrones = len(introns_sim)


		dict_anchors_finded = dict(list_anchors_finded)
		dict_anchors_sim = dict(list_anchors_sim)

		for i in set(introns_finded):

			ichr = i.split(':')[0]		
			istart = 0
			iend = 0
			strand = ''
		
			if '+' in i:
				istart = int(i.split(':')[1].split('+')[0])
				iend = int(i.split(':')[1].split('+')[1])
				strand = '+'
			elif '-' in i:
				istart = int(i.split(':')[1].split('-')[0])
				iend = int(i.split(':')[1].split('-')[1])
				strand = '-'

			dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
			if strand == '-':
				dn = dn.reverse_complement()

			#print row[3], i, str(dn).upper(), DRU_list, DRD_list, intron1_DR, intron2_DR, intron3_DR

			if i in OK_finded and i!='':
				N_intrones -= 1
				anchor = dict_anchors_finded[i]
				print row[3], i, read_length ,anchor, "OK", str(dn).upper()   #, DRU_list, DRD_list, intron1_DR, intron2_DR, intron3_DR

			elif i in Wrong and i!='':
				N_intrones -= 1
				anchor = dict_anchors_finded[i]
				print row[3], i, read_length , anchor,  "Wrong", str(dn).upper()  #, DRU_list, DRD_list, intron1_DR, intron2_DR, intron3_DR

		if N_intrones > 0:
			for i in introns_sim:

				ichr = i.split(':')[0]		
				istart = 0
				iend = 0
				strand = ''
		
				if '+' in i:
					istart = int(i.split(':')[1].split('+')[0])
					iend = int(i.split(':')[1].split('+')[1])
					strand = '+'
				elif '-' in i:
					istart = int(i.split(':')[1].split('-')[0])
					iend = int(i.split(':')[1].split('-')[1])
					strand = '-'

				dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]
				if strand == '-':
					dn = dn.reverse_complement()
			
				if not i in OK and N_intrones > 0 and i!='':
					N_intrones -= 1
					anchor = dict_anchors_sim[i]				
					print row[3], i, read_length , anchor, "NO_SJ_found", str(dn).upper()   #, DRU_list, DRD_list



    
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])            
	SJcheck(sys.argv[2])
