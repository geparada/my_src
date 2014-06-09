import sys
import csv
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq

'''Uso: python SJ_hist.py SJ.fasta SJ.sam SJ.check '''



def hist_anchor (introns, fasta, SAM, SJ_check):

	#Genome = dict(SeqTable)

	print >> sys.stderr, "Extrayendo informacion de intrones provenientes del trascriptoma modelo ..."

	reader3 = csv.reader(open(introns), delimiter = ' ')
	
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


		

	TOTAL_canonical = 0
	TOTAL_non_canonical = 0

	ambiguous_canonical = 0
	ambiguous_non_canonical = 0

	OK_canonical = 0
	OK_non_canonical = 0

	missed_canonical = 0
	missed_non_canonical = 0

	nosjfound_canonical = 0
	nosjfound_non_canonical = 0




	print >> sys.stderr, "Contando reads totales desde FASTQ ..."

        f = open(fasta)

        for record in SeqIO.parse(f, "fastq"):
                read = str(record.id)	
		introns_totales_FASTQ = read.split('=')[1].split('<>')
		for i in introns_totales_FASTQ:
			try:
				if type_intron[i] == 'Canonical':
					TOTAL_canonical += 1
				else:
					TOTAL_non_canonical += 1
			except KeyError:     
				pass

				#Los intrones de menos de 40pb no estan en la lista de intrones de genecode y no se consideraran para este analisis
	

	reader1 = csv.reader(open(SJ_check), delimiter = ' ')
	reader2 = csv.reader(open(SAM), delimiter = '\t')


	print >> sys.stderr, "Contantdo alineamientos ambiguos desdel el SAM..."

	for row in reader2:
		try:
			if int(row[13].split(':')[2]) == 2:
				introns_totales_SAM = row[0].split('=')[1].split('<>')
				for i in introns_totales_SAM:
					try:
						if type_intron[i] == 'Canonical':
							ambiguous_canonical += 1
						else:
							ambiguous_non_canonical += 1
					except KeyError:     
						pass

		except IndexError:             #Permite saltarse el header del SAM
			pass
			
	print >> sys.stderr, "Contantdo por tipo de error desde SJ.check..."

	for row in reader1:
		i = row[1]
		status = row[4]
		error_type = ",".join(row[5:8]).replace(",", " ")
		

		if status == 'OK':
			try:
				if type_intron[i] == 'Canonical':
					OK_canonical += 1
				else:
					OK_non_canonical += 1
			except KeyError:     
				pass

		elif error_type != 'No SJ found':
			try:
				if type_intron[i] == 'Canonical':
					missed_canonical += 1
				else:
					missed_non_canonical += 1
			except KeyError:     
				pass
		elif error_type == 'No SJ found': 			
			try:
				if type_intron[i] == 'Canonical':
					nosjfound_canonical += 1
				else:
					nosjfound_non_canonical += 1
			except KeyError:     
				pass


	unmapped_canonical = (TOTAL_canonical - (OK_canonical + missed_canonical + nosjfound_canonical + ambiguous_canonical))
	unmapped_non_canonical = (TOTAL_non_canonical - (OK_non_canonical + missed_non_canonical + nosjfound_non_canonical + ambiguous_non_canonical))



	print 'Intron Type | Total | OK | Wrong | No SJ found | Ambiguous | Unmapped | %FP' 
				
	print 'Canonical', TOTAL_canonical, percent(OK_canonical, TOTAL_canonical) , percent(missed_canonical, TOTAL_canonical), percent(nosjfound_canonical, TOTAL_canonical), percent(ambiguous_canonical, TOTAL_canonical), percent(unmapped_canonical, TOTAL_canonical), percent(missed_canonical, missed_canonical + OK_canonical)

	print 'Non-Canonical', TOTAL_non_canonical, percent(OK_non_canonical, TOTAL_non_canonical) , percent(missed_non_canonical, TOTAL_non_canonical), percent(nosjfound_non_canonical, TOTAL_non_canonical), percent(ambiguous_non_canonical, TOTAL_non_canonical), percent(unmapped_non_canonical, TOTAL_non_canonical), percent(missed_non_canonical, missed_non_canonical + OK_non_canonical)  				
  				


        f.close()
					

def percent (c, total):
	return (100*float(c))/float(total)






if __name__ == '__main__':
	#GenomicTabulator()
	hist_anchor (sys.argv[1], sys.argv[2],sys.argv[3],sys.argv[4])
		
		
