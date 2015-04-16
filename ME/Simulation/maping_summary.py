import sys
import csv
from collections import defaultdict 

total_count = defaultdict(int)

total_reads = 2458254833

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def main(log):
	for row in csv.reader(open(log), delimiter = ':'):
		if len(row)==2:
			total_count[row[0]] += int(row[1])




	Number_of_mapped_reads = total_count['Number of mapped reads']

	Number_of_unmapped_reads = total_count['Number of unmapped reads']

	Number_of_reads_after_preprocessing = Number_of_mapped_reads + Number_of_unmapped_reads

	Number_of_uniq_mapped_reads = total_count['Number of uniq mapped reads']	

	Number_of_uniq_SJ_mapped_reads = total_count['Number of uniq SJ mapped reads']






	print "Total reads", total_reads, percent(total_reads, total_reads)
	print "Number of reads after preprocessing", Number_of_reads_after_preprocessing, percent(Number_of_reads_after_preprocessing, total_reads)
	print "Number of mapped reads", Number_of_mapped_reads, percent(Number_of_mapped_reads, total_reads)
	print 'Number of uniq mapped reads', Number_of_uniq_mapped_reads, percent(Number_of_uniq_mapped_reads, total_reads)
	print 'Number of uniq SJ mapped reads', Number_of_uniq_SJ_mapped_reads, percent(Number_of_uniq_SJ_mapped_reads, total_reads)


if __name__ == '__main__':
	main(sys.argv[1])