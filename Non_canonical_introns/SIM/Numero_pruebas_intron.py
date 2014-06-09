import sys
import csv
from collections import defaultdict
from operator import itemgetter


def percent (c, total):
        return (100*float(c))/float(total)


def Intron_proof(bed12):
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
    
	count_intrones_OK = defaultdict(int)     #contar distinto
	count_intrones_missed = defaultdict(int)
	count_intrones_existent = defaultdict(int)
	count_intrones_nonexistent = defaultdict(int)

	frec_OK = defaultdict(int)
	frec_missed = defaultdict(int)
	frec_existent = defaultdict(int)
	frec_nonexistent = defaultdict(int)	


	keys_OK = []
	keys_missed = []
	keys_existent = []
	keys_nonexistent = []
    
	introns_finded_total = []
	introns_sim_total = []
   
	for row in reader1:
      
		q = map (int, row[11].strip(",").split(","))                      # obteniendo lista int de las distintas variables de la linea psl
		b = map(int, row[10].strip(",").split(","))
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
		introns_finded_total += introns_finded

		print introns_finded_total
	        
		introns_sim = row[3].split('=')[1].split('<>')
		#introns_sim_total +=  introns_sim                             
		
	
		for i in introns_finded:
			if i in introns_sim and i!='':
				count_intrones_OK[i] += 1
				keys_OK.append(i)
	
			elif i!='':
				count_intrones_missed[i] += 1
				keys_missed.append(i)




	for i in introns_finded_total:
		if count_intrones_OK.has_key(i) == True and i!='':
			count_intrones_existent[i] += 1
			keys_existent.append(i)
	
		elif i!='':
			count_intrones_nonexistent[i] += 1
			keys_nonexistent.append(i)

		

	

	for key in set(keys_OK):
		frec_OK[ count_intrones_OK[key] ] += 1

        for key in set(keys_missed):
                frec_missed[ count_intrones_missed[key] ] += 1

        for key in set(keys_existent):
                frec_existent[ count_intrones_existent[key] ] += 1

        for key in set(keys_nonexistent):
                frec_nonexistent[ count_intrones_nonexistent[key] ] += 1

	print "Proof_number", "OK", "Wrong", "existent", "nonexistent", "%FP"

	for n in range(1,101):
		
		OK = frec_OK[n]
		Wrong = frec_missed[n]
		Existent = frec_existent[n]
		Nonexistent = frec_nonexistent[n]
			
		try: 
			print n, OK, Wrong, Existent, Nonexistent, percent(Wrong, OK + Wrong)

		except ZeroDivisionError:
			
			print n, OK, Wrong, Existent, Nonexistent, 0 

#        print count_intrones_missed

		
			

    
if __name__ == '__main__':            
	 Intron_proof(sys.argv[1])
