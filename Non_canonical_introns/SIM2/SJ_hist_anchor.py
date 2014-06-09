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



def hist_anchor (fasta, SAM, SJ_check):

	hist_OK = defaultdict(int)
	hist_missed = defaultdict(int)
	hist_nosjfound = defaultdict(int)
	hist_ambiguous = defaultdict(int)
	hist_TOTAL = defaultdict(int)


	#sorted_OK = sorted(hist_OK.iteritems(), key=itemgetter(1))
	#sorted_missed = sorted(hist_missed.iteritems(), key=itemgetter(1))
	#sorted_nosjfound = sorted(hist_nosjfound.iteritems(), key=itemgetter(1))

	print >> sys.stderr, "Contando reads totales desde FASTQ ..."

        f = open(fasta)

        for record in SeqIO.parse(f, "fastq"):
                read = str(record.id)	
		anchors_totales_FASTQ = read.split('=')[2].split('<>')
		for a in anchors_totales_FASTQ:
			hist_TOTAL[min(map(int, a.split('-')))] += 1 





	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	reader2 = csv.reader(open(SAM), delimiter = '\t')


	print >> sys.stderr, "Contantdo alineamientos ambiguos desdel el SAM..."

	for row in reader2:
		try:
			if int(row[13].split(':')[2]) == 2:
				anchors_totales_SAM = row[0].split('=')[2].split('<>')
				for a in anchors_totales_SAM:
					hist_ambiguous[min(map(int, a.split('-')))] += 1


		except IndexError:
			pass
			
	print >> sys.stderr, "Contantdo por tipo de error desde SJ.check..."

	for row in reader1:
		status = row[4]
		error_type = ",".join(row[5:8]).replace(",", " ")
		

		if status == 'OK':
			hist_OK[row[3]] += 1

		elif error_type != 'No SJ found':
			hist_missed[row[3]] += 1

		elif error_type == 'No SJ found': 			
			hist_nosjfound[row[3]] += 1

	print 'Anchor | Total | OK | Wrong | No SJ found | Ambiguous | Unmapped | %FP' 

	for nt in range(1,51):                   #column join
		OK = hist_OK[str(nt)]
		Wrong = hist_missed[str(nt)]
		NO_SJ = hist_nosjfound[str(nt)]
		Total = hist_TOTAL[nt]
		Ambiguous = hist_ambiguous[nt]
		Unmapped = (Total - (OK + Wrong + NO_SJ + Ambiguous))

				
		print nt, Total, percent(OK, Total) , percent(Wrong, Total), percent(NO_SJ, Total), percent(Ambiguous, Total), percent(Unmapped, Total), percent(Wrong, Wrong + OK)  				


        f.close()
					

def percent (c, total):
	return (100*float(c))/float(total)






if __name__ == '__main__':
	#ReadsTabulator(sys.argv[1])
	hist_anchor (sys.argv[1],sys.argv[2],sys.argv[3])
		
		
