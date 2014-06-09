import sys
import csv
import random 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable=[]

def Genomictabulator(fasta):
	
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
        
		q = map (int, row[11].strip(",").split(","))                      # obteniendo lista int de las distintas variables de la linea psl
		b = map(int, row[10].strip(",").split(","))
		start = int(row[1])
		bn = int(row[9])	                  #El blocknum es importante para saber cuantos intrones hay en el alineamiento
			                                                   
		intron1 = ''
		intron2 = ''
		intron3 = ''
	
	

		if bn>=2 and ((start + q[1]) - (start + q[0] + b[0]))>=40:
			intron1 = row[0] + ':' +  str(start + q[0] + b[0]) + row[5] + str(start + q[1])
   
		if bn>=3 and ((start + q[2]) - (start + q[1] + b[1]))>=40:
			intron2 = row[0] + ':' +  str(start + q[1] + b[1]) + row[5] + str(start + q[2])

		if bn>=4 and ((start + q[3]) - (start + q[2] + b[2]))>=40:
			intron3 = row[0] + ':' +  str(start + q[2] + b[2]) + row[5] + str(start + q[3])
        
		introns_finded = [intron1, intron2, intron3]
		introns_finded_DR = [intron1, intron2, intron3]	        

		introns_sim = row[3].split('=')[1].split('<>')
		anchors_sim = row[3].split('=')[2].split('<>')                                   
	
		bn_sim = len(introns_sim) + 1

		read_length = sum(map(int, anchors_sim[0].split('-')))
	
		if bn_sim==3:
			read_length = read_length + int(anchors_sim[1].split('-')[1])
		if bn_sim==4:
			read_length = read_length + int(anchors_sim[1].split('-')[1]) + int(anchors_sim[2].split('-')[1]) 	

#		print row


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

#			print i, istart, iend
		#	introns_DR = [i]

			for chr in SeqTable:
				L = 100         #Solo permite que se corra L pares de bases para buscar DR 		

				if chr[0] == ichr:
					SJ5U = chr[1][istart-L : istart].lower()
					SJ5D = chr[1][istart : istart+L].lower()
					SJ3U = chr[1][iend-L : iend].lower()
					SJ3D = chr[1][iend : iend+L].lower()
				
					if '-' in i:
						SJ5U = chr[1][iend : iend+L].lower().reverse_complement() 
						SJ5D = chr[1][iend-L : iend].lower().reverse_complement()
						SJ3U = chr[1][istart : istart+L].lower().reverse_complement()
						SJ3D = chr[1][istart-L : istart].lower().reverse_complement()
				
					DRU = 0
					DRD = 0
										
				
					try:
						while SJ5U[L-1-DRU]==SJ3U[L-1-DRU]:
							DRU = (DRU + 1)
							introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart-DRU) + strand + str(iend-DRU)]
							if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
								break
					except IndexError:
						pass 

					try:
						while SJ5D[DRD]==SJ3D[DRD]:
							DRD = (DRD + 1)
							introns_finded_DR = introns_finded_DR + [ichr + ':' + str(istart+DRD) + strand + str(iend+DRD)]
							if SJ5D[DRD]!=SJ3D[DRD]:
								break
					except IndexError:
						pass

#					print row[3], SJ5U[-20:], SJ5D[:20], SJ3U[-20:], SJ3D[:20], DRU, DRD, introns_finded_DR

						

	
	
	
		for i, a in zip(introns_sim, anchors_sim):
					
			#print row[3], introns_finded_DR

			if i in introns_finded_DR:
				print row[3], i, read_length, min(map(int, a.split('-'))), "OK", ",".join(introns_finded_DR)
			
			elif intron1=='' :
				print row[3], i, read_length, min(map(int, a.split('-'))), "ERROR", "No SJ found"
	
			elif intron1!='' :
				print row[3], i, read_length, min(map(int, a.split('-'))), "ERROR", ",".join(introns_finded_DR)






    
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])            
	SJcheck(sys.argv[2])
