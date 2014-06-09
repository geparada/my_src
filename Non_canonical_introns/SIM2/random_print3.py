import sys
import random
import csv

def ran_printer(bigfile, total_lines, lines_wanted):

	reader = csv.reader(open(bigfile), dialect='excel-tab' )
	
	removed = 0
	current_len = total_lines - removed
	to_remove = total_lines - lines_wanted
	

	while to_remove > 0:
		for row in reader:
			probability = float(to_remove) / float(current_len )
			if random.random() < probability:
            			to_remove -= 1
				current_len -= 1
			else:
				current_len -= 1
				print " ".join(row)
			

if __name__ == '__main__':
	ran_printer(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
