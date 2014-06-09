import sys
import csv


def main(sam,bed12):
	read_seq_list = []

	reader1 = csv.reader(open(sam), delimiter = '\t')	
	reader2 = csv.reader(open(bed12), delimiter = '\t')

	
	print >> sys.stderr, "Extranyendo secuencia de los reads desde el SAM ...",
	for row in reader1:
		
		
		if row[0] != '@SQ':
			read = row[0]
			CIGAR = row[5]
			seq = row[9]

			if 'N' in CIGAR: 
				read_seq_list.append((read, seq))
	
	read_seq_dict = dict(read_seq_list)   #diccionario con el nombre del read y la secuencia
	
	print >> sys.stderr, "OK"
	
	for row in reader2:
		read = row[3]
		seq = read_seq_dict[read].upper()
		bedseq = row + [seq]
		
		print '\t'.join(bedseq)
	
	


if __name__ == '__main__':
	
	main(sys.argv[1],sys.argv[2]) 
