import sys
import csv
import random

RNA = []
Genoma = []
SeqName = []
Resultados=[]


def SJtag(AliSeq, range1, range2):
    '''Genera reads con splice juntios desde un transcriptoma de referencia

    *Input file = USCS BLAT alignments table (21 columns ) joined with fasta converted to tab (3 columns). This file has 24 columns.
    *Usage = python file longitud_minima longitud_maxima numero_de_ciclos
    *Output = seq ID , intron name, strand, tag sequence, tag length '''       
        

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
        
        qstarts = map (int, row[19].split(",")[1:-1])                      # obteniendo lista int de las distintas variables de la linea psl
        tstarts = map(int, row[20].split(",")[0:-1])
        blocksizes = map(int, row[18].split(",")[0:-1])
       
	
	for t0, b0, t1, b1, t2, b2, t3, b3, q in zip(tstarts, blocksizes, tstarts[1:], blocksizes[1:], tstarts[2:], blocksizes[2:], tstarts[3:], blocksizes[3:], qstarts[1:]):

        	n1 = random.randrange(1,100)     #generando numeros alazar entre 1 y 100 para sumar/restar a los qstart
        	n2 = random.randrange(1,100)
	
            
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
         
                if range1<=len(tag) and range2>=len(tag) and (b0+b1)>=n1 and (b2+b3)>=n2:
                    #print row[13] + str(tbsum) + row[8] + str(t2) + row[9], row[13] + str(tbsum) + row[8] + str(t2) + tag
		    #print row[13] + str(tbsum) + row[8] + str(t2), tag, len(tag)
		    #print ">" + row[9] + "=" + row[13] + ":" + str(tbsum) + row[8] + str(t2) + "=" + str(n1) + "-" + str(n2) + "\n" + tag               
		    #print row[13] + str(tbsum) + row[8] + str(t2) + "\t" + row[9]
		    print ">" + row[9] + "=" + intron0 + intron1 + intron2 + "=" + anchor0 + anchor1 + anchor2 + "\n" + tag 
		    #print ">" + row[9] + "=" intron1 "="  anchor1 "\n" + tag	




def cycler(a, b, c, n):
        for i in range(n):
                SJtag(a, b, c)

    
if __name__ == '__main__':                                      #Permite tomar argumentos mediante el terminal
    cycler(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) )
