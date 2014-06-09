import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq

'''Uso: python SJ_hist.py SJ.fasta SJ.sam SJ.check '''

SeqTable=[]

def ReadsTabulator(fasta):

        f = open(fasta)

        for chrfa in SeqIO.parse(f, "fasta"):
                table = str(chrfa.id)   #, chrfa.seq
                SeqTable.append(table)

        f.close()



def hist_anchor ( SAM, SJ_check):
	hist_OK = defaultdict(int)
	hist_missed = defaultdict(int)
	hist_nosjfound = defaultdict(int)
	hist_unmapped = defaultdict(int)
	hist_TOTAL = defaultdict(int)


	#sorted_OK = sorted(hist_OK.iteritems(), key=itemgetter(1))
	#sorted_missed = sorted(hist_missed.iteritems(), key=itemgetter(1))
	#sorted_nosjfound = sorted(hist_nosjfound.iteritems(), key=itemgetter(1))

	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	reader2 = csv.reader(open(SAM), delimiter = '\t')

	for row in reader2:
		try:
			if row[2] == '*':
				anchors_totales = row[0].split('=')[2].split('<>')
				for a in anchors_totales:
					hist_unmapped[min(map(int, a.split('-')))] += 1


		except IndexError:
			pass
			

	for read in SeqTable:
	
		anchors_totales = read.split('=')[2].split('<>')
		for a in anchors_totales:
			hist_TOTAL[min(map(int, a.split('-')))] += 1 







	for row in reader1:
		status = row[4]
		error_type = ",".join(row[5:8]).replace(",", " ")
		

		if status == 'OK':
			hist_OK[row[3]] += 1

		elif error_type != 'No SJ found':
			hist_missed[row[3]] += 1

		elif error_type == 'No SJ found': 			
			hist_nosjfound[row[3]] += 1

	print 'Anchor','OK','Wrong','No SJ found', 'Ambiguous', 'Unmapped', "%FP" 

	for nt in range(1,51):                   #column join
		OK = hist_OK[str(nt)]
		Wrong = hist_missed[str(nt)]
		NO_SJ = hist_nosjfound[str(nt)]
		Unmapped = hist_unmapped[nt]
		Total = hist_TOTAL[nt]
		Ambiguous = (Total - (OK + Wrong + NO_SJ + Unmapped))
						
		print nt, percent(OK, Total) , percent(Wrong, Total), percent(NO_SJ, Total), percent(Ambiguous, Total), percent(Unmapped, Total), percent(Wrong, Wrong + OK)  				


					

def percent (c, total):
	return (100*float(c))/float(total)






if __name__ == '__main__':
	ReadsTabulator(sys.argv[1])
	hist_anchor (sys.argv[2],sys.argv[3])
		
		
