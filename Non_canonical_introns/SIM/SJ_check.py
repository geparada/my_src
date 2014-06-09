import sys
import csv




def SJcheck(bed12):
    '''Verifica si es que los reads simulados de SJ alinearon correctamente

    *Input file = SJ.bed12
    *Usage = python SJ_check file
    *Output =  '''       
        

    reader1 = csv.reader(open(bed12), dialect='excel-tab' )
         
         # row[0]  chr
         # row[1]  alignment start
         # row[2]  alignment end
         # row[3]  name (Refseq=intron0<>intron1=anchor0<>anchor1)
         # row[4]    
         # row[5]  strand
         # row[6]  aligment start
         # row[7]  aligment end
         # row[8]  
         # row[9] blocknum
         # row[10] blocksizes
         # row[11] qstarts         
    
    
   
    for row in reader1:
        
        q = map (int, row[11].split(",")[0:])                      # obteniendo lista int de las distintas variables de la linea psl
        b = map(int, row[10].split(",")[0:])
	start = int(row[1])
	bn = int(row[9])	
			                                                   
	intron1 = ''
	intron2 = ''
	intron3 = ''
	
	

	if bn>=2:
		intron1 = row[0] + ':' +  str(start + q[0] + b[0]) + row[5] + str(start + q[1])
   
	if bn>=3:
		intron2 = row[0] + ':' +  str(start + q[1] + b[1]) + row[5] + str(start + q[2])

	if bn>=4:
                intron3 = row[0] + ':' +  str(start + q[2] + b[2]) + row[5] + str(start + q[3])
        
	introns_finded = [intron1, intron2, intron3]
	        
	introns_sim = row[3].split('=')[1].split('<>')
	anchors_sim = row[3].split('=')[2].split('<>')                                   
	
	bn_sim = len(introns_sim) + 1

	read_length = sum(map(int, anchors_sim[0].split('-')))
	
	if bn_sim==3:
		read_length = read_length + int(anchors_sim[1].split('-')[1])
	if bn_sim==4:
		read_length = read_length + int(anchors_sim[1].split('-')[1]) + int(anchors_sim[2].split('-')[1]) 	

	
	
	for i, a in zip(introns_sim, anchors_sim):
		if i in introns_finded:
			print row[3], i, read_length, min(map(int, a.split('-'))), "OK", intron1, intron2, intron3

		elif intron1=='' :
			print row[3], i, read_length, min(map(int, a.split('-'))), "ERROR", "No SJ found"
	
		elif intron1!='' :
                        print row[3], i, read_length, min(map(int, a.split('-'))), "ERROR", intron1, intron2, intron3


    
if __name__ == '__main__':            
	 SJcheck(sys.argv[1])
