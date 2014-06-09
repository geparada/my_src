import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq

'''Uso: python SJ_hist.py SJ.fasta SJ.sam SJ.check '''

SeqTable=[]

#def ReadsTabulator(fasta):

#        f = open(fasta)

#        for chrfa in SeqIO.parse(f, "fastq"):
#                table = str(chrfa.id)   #, chrfa.seq
#                SeqTable.append(table)

#        f.close()



def hist_anchor (fastq_stats, SAM, SJ_check):

	hist_OK = defaultdict(int)
	hist_missed = defaultdict(int)
	hist_nosjfound = defaultdict(int)
	hist_ambiguous = defaultdict(int)
	hist_unmapped = defaultdict(int)


	#sorted_OK = sorted(hist_OK.iteritems(), key=itemgetter(1))
	#sorted_missed = sorted(hist_missed.iteritems(), key=itemgetter(1))
	#sorted_nosjfound = sorted(hist_nosjfound.iteritems(), key=itemgetter(1))

	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	reader2 = csv.reader(open(SAM), delimiter = '\t')
	reader3 = csv.reader(open(fastq_stats), delimiter = ' ')

	stats = []

	for row in reader3:
		stats.append((row[0], row[2]))   #row[1] para canonicos y row[2] para no canonicos
	
	fastq = dict(stats)		


	print >> sys.stderr, "Contando alineamientos no mapeados y ambiguos desdel el SAM..."

	for row in reader2:
		try:
			anchors_totales_SAM = row[0].split('=')[2].split('<>')
			for a in anchors_totales_SAM:
				if row[2] == '*':
					hist_unmapped[min(map(int, a.split('-')))] += 1
				elif int(row[13].split(':')[2]) == 2:
					hist_ambiguous[min(map(int, a.split('-')))] += 1



		except IndexError:
			pass
			
	print >> sys.stderr, "Contando por tipo de error desde SJ.check..."

	for row in reader1:
		status = row[4]

		if status == 'OK':
			hist_OK[row[3]] += 1

		elif status == 'Wrong':
			hist_missed[row[3]] += 1

		elif status == 'NO_SJ_found': 			
			hist_nosjfound[row[3]] += 1

	print 'Anchor | Total FASTQ | Total SAM | OK | Wrong | No SJ found | Ambiguous | Unmapped | %True discovery rate | %Sensitivity'

	for nt in range(1,51):                   #column join
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
		



				
		print nt, Total_fastq, Total_sam, percent(OK, Total) , percent(Wrong, Total), percent(NO_SJ, Total), percent(Ambiguous, Total), percent(Unmapped, Total), percent(OK, Wrong + OK), percent(OK, Total_fastq)		


					

def percent (c, total):
	return (100*float(c))/float(total)






if __name__ == '__main__':
	hist_anchor (sys.argv[1],sys.argv[2],sys.argv[3])
		
		
