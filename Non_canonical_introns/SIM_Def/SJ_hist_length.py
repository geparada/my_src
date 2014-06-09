import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq

'''Uso: python SJ_hist.py SJ.fasta SJ.sam SJ.check '''


def hist_length(fastq_stats, SAM, SJ_check):

	hist_OK = defaultdict(int)
	hist_missed = defaultdict(int)
	hist_nosjfound = defaultdict(int)
	hist_ambiguous = defaultdict(int)
	hist_unmapped = defaultdict(int)

	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	reader2 = csv.reader(open(SAM), delimiter = '\t')
	reader3 = csv.reader(open(fastq_stats), delimiter = ' ')

	stats = []

	for row in reader3:
		stats.append((row[0], row[1]))
	
	fastq = dict(stats)		


	print >> sys.stderr, "Contantdo alineamientos ambiguos desdel el SAM..."

	for row in reader2:


		try:

			anchors_sim = row[0].split('=')[2].split('<>')
			read_length = sum(map(int, anchors_sim[0].split('-')))
	
			if len(anchors_sim)==2:
				read_length = read_length + int(anchors_sim[1].split('-')[1])
			if len(anchors_sim)==3:
				read_length = read_length + int(anchors_sim[1].split('-')[1]) + int(anchors_sim[2].split('-')[1])
 
			introns_totales_SAM = row[0].split('=')[1].split('<>')
			for i in introns_totales_SAM:
				if row[2] == '*':
					hist_unmapped[read_length] += 1
				elif int(row[13].split(':')[2]) == 2:
					hist_ambiguous[read_length] += 1

		except IndexError:
			pass
			
	print >> sys.stderr, "Contantdo por tipo de error desde SJ.check..."

	for row in reader1:
		status = row[4]

		if status == 'OK':
			hist_OK[row[2]] += 1

		elif status == 'Wrong':
			hist_missed[row[2]] += 1

		elif status == 'NO_SJ_found': 			
			hist_nosjfound[row[2]] += 1

	print 'Length | Total FASTQ | Total SAM | OK | Wrong | No SJ found | Ambiguous | Unmapped | %True discovery rate | %Sensitivity'

	for nt in range(1,101):                   #column join
		OK = hist_OK[str(nt)]
		Wrong = hist_missed[str(nt)]
		NO_SJ = hist_nosjfound[str(nt)]
		Ambiguous = hist_ambiguous[nt]
		Unmapped = hist_unmapped[nt]
		Total_sam = OK + Wrong + NO_SJ + Ambiguous + Unmapped
		Total_fastq = fastq[str(nt)]
		Total = Total_fastq

		if int(Total_fastq) >= int(Total_sam):
			
			Unmapped = int(Unmapped) + (int(Total_fastq) - int(Total_sam))

		else:
			Total = int(Total_sam)
				
		if Total!=0 and Total!='':
			print nt, Total_fastq, Total_sam, percent(OK, Total) , percent(Wrong, Total), percent(NO_SJ, Total), percent(Ambiguous, Total), percent(Unmapped, Total), percent(OK, Wrong + OK), percent(OK, Total_fastq)	
		else:
			print nt, 0, 0, 0, 0, 0, 0, 0 				


def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0





if __name__ == '__main__':
	#ReadsTabulator(sys.argv[1])
	hist_length(sys.argv[1],sys.argv[2],sys.argv[3])
		
		
