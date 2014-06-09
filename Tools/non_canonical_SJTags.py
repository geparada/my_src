import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

SeqTable1 = []
SeqTable2 = []
gene_model = defaultdict(list)

#USE: python ~/my_src/Tools/non_canonical_SJTags.py $Genoma $Model_fasta $Model_bed12 $pre_final_table


def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		table = str(chrfa.id), chrfa.seq
		SeqTable1.append(table)

	f.close()  

	print >> sys.stderr, "OK"

def Transcriptometabulator(refseq_fasta):
	
	print >> sys.stderr, "Cargando transcriptoma en la memoria RAM ...",	
	
	for record in SeqIO.parse(refseq_fasta, "fasta"):
		ID = record.id.split('|')[0]
		table = (str(ID), record.seq)
		SeqTable2.append(table)

	print >> sys.stderr, "OK"

def gene_model_extractor(bed12):

	reader = csv.reader(open(bed12), delimiter = '\t')
	
	for row in reader:
		chr = row[0]
		start = row[1]
		end = row[2]
		name = row[3]
		strand = row[5]
		block_sizes = row[10]
		tstarts = row[11]
		
		gene_model[chr + strand].append((name, start, end, strand, block_sizes, tstarts))
		
		

def main(pre_final_table):
	Genome = dict(SeqTable1)
	Transcriptome = dict(SeqTable2)
	
#	intron_list = set([])

	csv.field_size_limit(1000000000)
	
	reader = csv.reader(open(pre_final_table), delimiter = '\t')
	
	
	for row in reader:
		ichr = row[0]
		istart = int(row[1])
		iend = int(row[2])
		istrand = row[5]
		intron = ichr + ':' + str(istart) + istrand + str(iend)
		dn = row[3].split("|")[1]


		
		
		if dn!='GTAG' and dn!='GCAG' and dn!='ATAC':
			
			tag_genomic_up = Genome[ichr][istart-100:istart]
			tag_genomic_down = Genome[ichr][iend:iend+100]		
					

			for transcript in gene_model[ichr + istrand]:
				qname = transcript[0]
				start = int(transcript[1])
				end = int(transcript[2])
				strand = transcript[3]
				block_sizes = map (int, transcript[4].strip(",").split(","))		
				tstarts_0 = map (int, transcript[5].strip(",").split(","))            #los llamo tstarts_0, por que estan referenciados a 0, para obtener los verdaderos tstarts hay que sumarle el start

				
								
				if start < istart and end > iend:
					
					
					SJ_up = []
					SJ_down = []
					
					qstart = 0
					qstarts = [0]                        #Es necesario para no generar errores de slice
					
					for t , b in zip(tstarts_0, block_sizes):
						qstart += b
						qstarts.append(qstart)
						
						
						up = start + t - istart
						down = start + t + b - iend
						
						if up <= 0:
							SJ_up.append(up)
						if down >= 0:
							SJ_down.append(down)
							
					
					genomic_slice_up = 0            #Coordenadas para hacerle un slice al tag de genomico si es que pasa por alguna SJ
					genomic_slice_down = 0
					
					if len(SJ_up) != 0:
						genomic_slice_up = -int(SJ_up[-1])        
					if len(SJ_down) != 0:
						genomic_slice_down = int(SJ_down[0])
										
					try:
					
						qstart_up = qstarts[len(SJ_up)-1]                 #Selecionando el bloque dentro del query para extraer secuencia desde el trasncrito
						qstart_down = qstarts[-len(SJ_down)]
					
						tag_complement_up = ''
						tag_complement_down = ''				
					
						tag_core_up = tag_genomic_up
						if genomic_slice_up < 100 and genomic_slice_up !=0 :
							tag_core_up = Genome[ichr][istart-genomic_slice_up:istart]
							tag_complement_up = Transcriptome[qname][qstart_up - (100 - genomic_slice_up) :qstart_up]
	#						if qstart_up - (100 - genomic_slice_up) > 0:
	#							tag_complement_up = Transcriptome[qname][:qstart_up]
							if strand == '-':
								tag_complement_up = Transcriptome[qname].reverse_complement()[qstart_up - (100 - genomic_slice_up) :qstart_up]
						
											
						tag_core_down = tag_genomic_down
						if genomic_slice_down < 100 and genomic_slice_down !=0:
							tag_core_down = Genome[ichr][iend:iend+genomic_slice_down]
							tag_complement_down = Transcriptome[qname][qstart_down: qstart_down + (100 - genomic_slice_down) ]
							if qstart_down + (100 - genomic_slice_down) > len(Transcriptome[qname]):
								tag_complement_down = Transcriptome[qname][qstart_down:]
							
							if strand == '-':
								tag_complement_down = Transcriptome[qname].reverse_complement()[qstart_down: qstart_down + (100 - genomic_slice_down) ]
						
						tag = 	(str(tag_complement_up) + str(tag_core_up) + str(tag_core_down) + str(tag_complement_down)).upper()
						if strand == '-':
							tag = str(Seq(tag).reverse_complement())
							TOTAL_block_up = str(len(str(tag_core_down) + str(tag_complement_down)))
							TOTAL_block_down = str(len(str(tag_complement_up) + str(tag_core_up)))
						
						ID = ">" + intron + "|" + qname + "|" + TOTAL_block_up + "_" + TOTAL_block_down	


						print ID
						print tag 
						
						
							
					except KeyError:
						pass
						
	
		
		
	
	


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	Transcriptometabulator(sys.argv[2])
	gene_model_extractor(sys.argv[3])
	main(sys.argv[4])
