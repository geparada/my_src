import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq

'''Uso: python SJ_hist.py SJ.fasta SJ.sam SJ.check '''



def hist_anchor (introns, fastq_stats, SAM, SJ_check):


	print >> sys.stderr, "Extrayendo informacion de intrones provenientes del trascriptoma modelo ..."

	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	reader2 = csv.reader(open(SAM), delimiter = '\t')
	reader3 = csv.reader(open(introns), delimiter = ' ')
	reader4 = csv.reader(open(fastq_stats), delimiter = ' ')
	
	Fastq_Canonical = 0
	Fastq_Non_canonical = 0
	for row in reader4:
		Fastq_Canonical = row[1]
		Fastq_Non_Canonical = row[2]

		
	
	list_info_introns = []

	for row in reader3:
		intron = row[1] + ':' + row[2] + row[4] + row[3]
		dn = row[7]
		dn_type = ''
		if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
			dn_type = 'Canonical'
		else:
			dn_type = 'Non-Canonical'
		list_info_introns.append((intron, dn_type))

	type_intron = dict(list_info_introns)

		

	unmapped_canonical = 0
	unmapped_non_canonical = 0

	ambiguous_canonical = 0
	ambiguous_non_canonical = 0

	OK_canonical = 0
	OK_non_canonical = 0

	missed_canonical = 0
	missed_non_canonical = 0

	nosjfound_canonical = 0
	nosjfound_non_canonical = 0


	print >> sys.stderr, "Contando alineamientos no mapeados y ambiguos desdel el SAM..."

	for row in reader2:
		try:
			introns_totales_SAM = row[0].split('=')[1].split('<>')
			for i in introns_totales_SAM:
				if row[2] == '*':
					try:
						if type_intron[i] == 'Canonical':
							unmapped_canonical += 1
						else:
							unmapped_non_canonical += 1
					except KeyError:     
						pass

				elif int(row[13].split(':')[2]) == 2:

					try:
						if type_intron[i] == 'Canonical':
							ambiguous_canonical += 1
						else:
							ambiguous_non_canonical += 1
					except KeyError:     
						pass
		except IndexError:
			pass 
	
	print >> sys.stderr, "Contantdo por tipo de error desde SJ.check..."

	for row in reader1:
		i = row[1]
		status = row[4]
		dn = row[5]

		

		if status == 'OK':
			try:
				if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
					OK_canonical += 1
				else:
					OK_non_canonical += 1
			except KeyError:     
				pass

		elif status == 'Wrong':            
			try:
				if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
					missed_canonical += 1
				else:
					missed_non_canonical += 1
			except KeyError:     
				pass

		elif status == 'NO_SJ_found': 			
			try:
				if dn == 'GTAG' or dn == 'GCAG' or dn == 'ATAC':
					nosjfound_canonical += 1
				else:
					nosjfound_non_canonical += 1
			except KeyError:     
				pass


	TOTAL_canonical = OK_canonical + missed_canonical + nosjfound_canonical + ambiguous_canonical + unmapped_canonical
	TOTAL_non_canonical = OK_non_canonical + missed_non_canonical + nosjfound_non_canonical + ambiguous_non_canonical + unmapped_non_canonical



	print 'Intron Type | Total FASTQ | Total SAM | OK | Wrong | No SJ found | Ambiguous | Unmapped | %True discovery rate | %Sensitivity' 
				
	print 'Canonical', Fastq_Canonical, TOTAL_canonical, percent(OK_canonical, TOTAL_canonical) , percent(missed_canonical, TOTAL_canonical), percent(nosjfound_canonical, TOTAL_canonical), percent(ambiguous_canonical, TOTAL_canonical), percent(unmapped_canonical, TOTAL_canonical), percent(OK_canonical, missed_canonical + OK_canonical), percent(OK_canonical ,Fastq_Canonical)

	print 'Non-Canonical', Fastq_Non_Canonical, TOTAL_non_canonical, percent(OK_non_canonical, TOTAL_non_canonical) , percent(missed_non_canonical, TOTAL_non_canonical), percent(nosjfound_non_canonical, TOTAL_non_canonical), percent(ambiguous_non_canonical, TOTAL_non_canonical), percent(unmapped_non_canonical, TOTAL_non_canonical), percent(OK_non_canonical, missed_non_canonical + OK_non_canonical), percent(OK_non_canonical ,Fastq_Non_Canonical)  				
  				

					

def percent (c, total):
	return (100*float(c))/float(total)


if __name__ == '__main__':
	hist_anchor (sys.argv[1], sys.argv[2],sys.argv[3],sys.argv[4])
		
		
