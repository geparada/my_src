import sys
import csv

def intron_extractor(bed12):

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

		
		introns_finded = [intron1, intron2, intron3]
	





if __name__ == '__main__':            
	 intron_extractor(sys.argv[1])
