import sys
import csv

def hms_to_minutes(t):
	h, m, s = [float(i) for i in t.split(':')]
	return 60*h + m + s/60
	
def main(file):
	reader = csv.reader(open(file), delimiter = '\t')
	total_min = 0
	for row in reader:
		minutes = hms_to_minutes(row[2])
		total_min += minutes
	print total_min   

if __name__ == '__main__':
	main(sys.argv[1]) 
