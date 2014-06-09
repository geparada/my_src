import sys
import csv


def summary (SJ_check):

	reader = csv.reader(open(SJ_check), delimiter = ' ')
	OK = 0
	ERROR_missed = 0
	ERROR_nosjfound = 0


	for row in reader:
		status = row[4]
		error_type = ",".join(row[5:8]).replace(",", " ")
		

		if status == 'OK':
			OK = OK + 1		

		elif error_type == 'No SJ found':
			ERROR_nosjfound = ERROR_nosjfound + 1

		elif error_type != 'No SJ found': 			
			ERROR_missed = ERROR_missed + 1

	print  'OK SJ = ' + str(OK), 'Wrong SJ  = ' + str(ERROR_missed), 'No SJ found = ' + str(ERROR_nosjfound)





if __name__ == '__main__':
	summary(sys.argv[1])
		
		
