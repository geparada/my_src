import sys
import csv
import random

RNA = []
Genoma = []
SeqName = []
Resultados=[]
RefSeq = []


def expression (AliSeq):

	reader1 = csv.reader(open(AliSeq), dialect='excel-tab' )
		 
		# row[0]  matches
		# row[1]  misMatches
		# row[2]  repMaches
		# row[3]  nCount
		# row[4]  qNumInsert    
		# row[5]  qBaseInsert
		# row[6]  tNumInsert
		# row[7]  tBaseInsert
		# row[8]  strand
		# row[9] qName
		# row[10] qSize
		# row[11] qStart
		# row[12] qEnd         
		# row[13] tName
		# row[14] tSize
		# row[15] tStart
		# row[16] tEnd
		# row[17] blockCount
		# row[18] blockSizes
		# row[19] qStarts
		# row[20] tStarts

		# row[21] qName
		# row[22] 1
		# row[23] sequence

	for row in reader1:
		expression_factor = random.gauss(0.5, 0.2)
		
		RefSeq.append((row, expression_factor))
	


def SJtag(error_index):
	'''Genera reads con splice juntios desde un transcriptoma de referencia

	*Input file = USCS BLAT alignments table (21 columns ) joined with fasta converted to tab (3 columns). This file has 24 columns.
	*Usage = python file longitud_minima longitud_maxima numero_de_ciclos
	*Output = seq ID , intron name, strand, tag sequence, tag length '''

	reader2 = csv.reader(open(error_index), delimiter=' ' )
	
	nt = []
	nt_count = []
	P = []	
	

	for row in reader2:
		nt.append(row[0])
		nt_count.append(row[1])
		P.append(row[2])
	
			

	next_nt_count = nt_count[1:] + [(0)]
	Total = sum(map(float, nt_count)) - sum(map(float, next_nt_count))
	
	readlens = []	  #Lista de tamanos de reads que tienen la misma distribucion que el total
			
	for base, basecount, p, nextbasecount in zip(nt, nt_count, P, next_nt_count):
		readlencount = float(basecount) - float(nextbasecount)
		
		for x in range(int(readlencount/100)):      # se parte en 100 para que no se genere una lista tan pesada y asi optimizar memoria ram
			readlens.append(int(base)+1)
	
       
    
   
	for line in RefSeq:
		row, E_factor = line
		
		
        
		qstarts = map (int, row[19].split(",")[1:-1])                      # obteniendo lista int de las distintas variables de la linea psl
		tstarts = map(int, row[20].split(",")[0:-1])
		blocksizes = map(int, row[18].split(",")[0:-1])
       
	
		for t0, b0, t1, b1, t2, b2, t3, b3, q  in zip(tstarts, blocksizes, tstarts[1:], blocksizes[1:], tstarts[2:], blocksizes[2:], tstarts[3:], blocksizes[3:], qstarts[1:] ):

			#n1 = random.randrange(1,100)     #generando numeros alazar entre 1 y 100 para sumar/restar a los qstart
			#n2 = random.randrange(1,100)
			l = random.sample(readlens, 1)[0]
			
			n1 = random.randint(1,l-1)
			n2 = l - n1		
			
            
			tbsum0 = t0+b0
			tbsum1 = t1+b1
			tbsum2 = t2+b2
			#tbsum3 no es necesario
			                                                   
			start = q-n1                                                           
			end = q+n2

			intron0 = ''                                                  
			intron1 = row[13] + ":" + str(tbsum1) + row[8] + str(t2)                  
			intron2 = ''    

			anchor0 = ''
			anchor1 = str(n1) + "-" + str(n2)
			anchor2 = ''
		

# los intrones 0 y 2 estaran presente solo si el start o el end se escapan del rango en el cual se encuentran los exones mas proximos al SJ


			if n1>b1 and n2<=b2:
				intron0 = row[13] + ":" + str(tbsum0) + row[8] + str(t1) + "<>"

				anchor0 = str(n1-b1) + "-" + str(b1) + "<>"
				anchor1 = str(b1) + "-" + str(n2)

			if n2>b2 and n1<=b1:
				intron2 = "<>" + row[13] + ":" + str(tbsum2) + row[8] + str(t3)

				anchor1 = str(n1) + "-" + str(b2)
				anchor2 = "<>" + str(b2) + "-" + str(n2-b2)

			if n1>b1 and n2>b2:
				intron0 = row[13] + ":" + str(tbsum0) + row[8] + str(t1) + "<>"
				intron2 = "<>" + row[13] + ":" + str(tbsum2) + row[8] + str(t3)

				anchor0 = str(n1-b1) + "-" + str(b1) + "<>"
				anchor1 = str(b1) + "-" + str(b2)
				anchor2 = "<>" + str(b2) + "-" + str(n2-b2)


                    
			if row[8] == '+' :                          #Para los que aliniean en la hebra +
    
                    
				if start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
					start = 0  
				if end>len(row[23]):
					end=len((row[23]))
				tag = row[23][start:end]
    
                    
			if row[8] == '-' :                                #Para los que alinian en la hebra - es todo al inverso
    
				if end>len(row[23]):
					end=len(row[23])                                       
				tag = row[23][-end:-start]
				if start<=0: 
					tag = row[23][-end:]
			
			E_chance = random.random()

         
			if E_factor >= E_chance and (b0+b1)>=n1 and (b2+b3)>=n2:              #asegura que solo tengan 3 SJ como maximo
				print ">" + row[9] + "=" + intron0 + intron1 + intron2 + "=" + anchor0 + anchor1 + anchor2 + "\n" + tag 




def cycler(a, n):
	for i in range(n):
		SJtag(a)

    
if __name__ == '__main__':
	expression(sys.argv[1])                        
	cycler(sys.argv[2], int(sys.argv[3]) )
