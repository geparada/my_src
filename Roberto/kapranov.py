import sys
import csv

def main (with_Ts_sam, without_Ts_sam):
	""" Compara los dos SAMs y deja los que no alinearon al genoma con las Ts, pero si sin las Ts """
	
	with_Ts_unmapped = set([])

	for row in csv.reader(open(with_Ts_sam), delimiter = '\t') :

		if row[0][0]!="@":  # Para saltarse el header
			
			read = row[0]
			flag = row[1] 
			
			if flag=="4":
				with_Ts_unmapped.add(read)
	
	C = 0
		
	for row in csv.reader(open(without_Ts_sam), delimiter = '\t') :
		
		if row[0][0]!="@":  # Para saltarse el header
		
			read = row[0]
			flag = row[1] 
			
			if (read in with_Ts_unmapped) and (flag=="0" or flag=="16"):
				C += 1
				print "\t".join(row)
		
		else:
			print "\t".join(row)
	
	print >> sys.stderr, "There are " + str(C) + " alinments that only map without Ts"
	
				
if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])  
