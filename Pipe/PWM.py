import sys
import csv

def main(headers, data):
	reader1 = csv.reader(open(headers), delimiter = '\t')
	reader2 = csv.reader(open(data), delimiter = '\t')
	
	
	IDs = []
	
	for row in reader1:
		ID = row[0]
		IDs.append(ID)

	counter = 0
	
	for row in reader2:
		pos = row[0]
		A = row[1]
		C = row[2]		
		G = row[3]		
		T = row[4]
		
		if pos == "Pos":
			ID = IDs[counter]
			counter += 1
			print ">" + ID		
			
		else:
			print "\t".join([A, C, G, T])
			
	







if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])
