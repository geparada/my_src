import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
import re

genomic = open("temporal_genomic.sam", 'w')

sense_SJ = open("sense_SJ.sam", 'w')
antisense_SJ = open("antisense_SJ.sam", 'w')
#sense_SJ_uniq = open("sense_SJ_uniq.sam", 'w')
#antisense_SJ_uniq = open("antisense_SJ_uniq.sam", 'w')

Genome = {}

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()

def percent (c, total):
        try:
                return (100*float(c))/float(total)
        except ZeroDivisionError:
                return 0


def main(small_RNAs_SAM, small_RNAs_SAM_adapter, tags_fasta, refseq_fasta, refseq_gene_names, reflink, refseq_bed):


	ref_chr = {}
	ref_gene = {}
	refseq_seq = {}
	tags_seq = {}
	not_CG = set([])
	genomic_adapter = set([])
	genomic_adapter_SJ = set([])
	
	for record in SeqIO.parse(tags_fasta, "fasta"):
		tags_seq[record.id] = str(record.seq).upper()

	for record in SeqIO.parse(refseq_fasta, "fasta"):
		refseq_seq[record.id] = str(record.seq).upper()
			
	for row in csv.reader(open(small_RNAs_SAM_adapter), delimiter = '\t'):
		
		read = row[0]
		flag = row[1]
		tname = row[2]
		start = row[3] 		
		MD_tag = row[12]
		
		if flag == "0":
			start = str(int(start) + 6)
		
		first_mistmatch = ""

		i = 0		
		if flag =="16":
			MD_tag = MD_tag[::-1]
			
			
		while i < len(MD_tag):
			if MD_tag[i] in "0123456789":
				first_mistmatch += MD_tag[i]
			if MD_tag[i] in "ACTGN":
				break
			i+=1
		
		if flag =="16":
			first_mistmatch = first_mistmatch[::-1]
			
		first_mistmatch = int(first_mistmatch)
		
		if first_mistmatch >= 6:
			genomic_adapter.add(read+tname+start)
			genomic_adapter_SJ.add(read+tname)

	for row in csv.reader(open(refseq_bed), delimiter = '\t'):
		chr = row[0]		
		ID = row[3]		
		ref_chr[ID] = chr			
	
	for row in csv.reader(open(refseq_gene_names), delimiter = '\t'):
		ID = row[0]
		chr = row[1]
		gene = row[2]
		ref_gene[ID] = gene
	
	for row in csv.reader(open(reflink), delimiter = '\t'):
		ID = row[2]
		gene = row[0]
		ref_gene[ID] = gene					


	previous_row = []
	previous_read = ""
	reads_rows = []

	
	for row in csv.reader(open(small_RNAs_SAM), delimiter = '\t') :

		if row[0][0]!="@":  # Para saltarse el header
			
			read = row[0]
			flag = row[1] 
			
			if flag=="0" or flag=="16":


				read = row[0]
				ali_flag = row[1]							
				start = int(row[3])           #Sam es 1 referenciado 
				MAPQ = row[4]
				cigar = row[5]
				seq = row[9]
				qual = row[10]
				tname = row[2]
				matches = int(cigar.strip('M'))															
				anchor_up =  92- start + 1
				anchor_down =  start + matches - 92 - 1
		
				mismatches = int(row[13].strip("NM:i:"))
				MD_tag = row[12]
				NM_tag = row[13]
						
													
				if "|" in tname: # alinea a un SJ
										
					try:
									
						ref_ID, intron = tname.split("|")
									
						#Cambiando hebra
									
						new_flag = "0"
						
						if (("+" in intron) and ali_flag=="16") or (("-" in intron) and ali_flag=="0"):
							new_flag = "16"																					
									
						#Coordenadas de intron
									
						chr = ref_chr[ref_ID]

						chrstart, end = re.findall(r"[\w']+", intron)  #no tiene : asi que se parte en 2 solamente
						start = ""
									
																		
						if (chr in chrstart):
							start = chrstart.replace(chr, "")
										
						elif (chr.split("_")[0] in chrstart):
							chr = chr.split("_")[0]
							start = chrstart.replace(chr, "")

						elif ("X" in chrstart):
							chr = "chrX"
							start = chrstart.replace(chr, "")

						elif ("Y" in chrstart):
							chr = "chrY"
							start = chrstart.replace(chr, "")
										
						elif (end[:3] in chrstart):
							chr = chrstart.split(end[:3])[0]
							start = chrstart.replace(chr, "")

						elif (end[:2] in chrstart):
							chr = chrstart.split(end[:2])[0]
							start = chrstart.replace(chr, "")
						
						elif (end[:1] in chrstart):
							chr = chrstart.split(end[:1])[0]
							start = chrstart.replace(chr, "")																			
										
	
						ilength = int(end) -int(start)
									
						new_cigar =  str(anchor_up) + "M" + str(ilength)+  "N" + str(anchor_down) + "M"											

						if anchor_up >= 6 and anchor_down >= 6:
										
							ali_start = int(start) - anchor_up + 1										

							XN_tag = "XN:Z:?"
							XI_tag = "XI:Z:" + intron
											
							if ref_gene.has_key(ref_ID):
								XN_tag = "XN:Z:" + ref_gene[ref_ID]

							XS_tag = "XS:A:+"
							new_MD_tag = MD_tag 

							new_seq = seq
							new_phred = qual	
							if	new_flag != ali_flag:
								new_seq = str(Seq(seq).reverse_complement())
								new_phred = qual[::-1]
								new_cigar = str(anchor_down) + "M" + str(ilength)+  "N" + str(anchor_up) + "M"
								ali_start = int(start) - anchor_down + 1
											
								MD_number = []
								MD_letter = ""									
#								for i in MD_tag[::-1]:
#									if i in set("ACGTN"):
#										MD_letter = i
#										new_MD_tag += "".join(MD_number[::-1]) + MD_letter
#										MD_number = []
#									elif i in {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}:
#										MD_number.append(i)
#									elif i == "Z":
#										new_MD_tag += "".join(MD_number[::-1])
										

							XP_tag = "XP:Z:F"	#poliA tag
							XM_tag = "XM:Z:F"	#Adaptador genomico											

											
							if read+tname in genomic_adapter_SJ:
								XM_tag = "XM:Z:T"											
										
							if (refseq_seq.has_key(ref_ID)) and (tags_seq.has_key(tname)):
								'''-Se splitea el refseq con la secuencia del tag para buscar poliAs 'AAAAAA 3''
								   -Si es que el refseq tiene mas de un segmento con la misma secuencia del tag, entonces se revisan todos los posibles segmentos 3' ''' 
								tag_up, tag_down = tags_seq[tname][:len(tags_seq[tname])/2], tags_seq[tname][len(tags_seq[tname])/2:]
											
								if (tags_seq[tname] in refseq_seq[ref_ID]):
												
												
									if ali_flag == "0":
											
										for end3 in refseq_seq[ref_ID].split(tags_seq[tname])[1:]:

											if matches < 95:
												if 'AAAAAA' in (tag_down + end3)[anchor_down:anchor_down+10]:
													XP_tag = "XP:Z:T"
											else:
												if 'AAAAAA' in (tag_down + end3)[anchor_down:anchor_down+120]:
													XP_tag = "XP:Z:T"

									if ali_flag == "16":
											
										for end3 in refseq_seq[ref_ID].split(tags_seq[tname])[:-1]:

											if matches < 95:
												if 'TTTTTT' in (end3 + tag_up)[-anchor_up-10:-anchor_up]:
													XP_tag = "XP:Z:T"
											else:
												if 'TTTTTT' in (end3 + tag_up)[-anchor_up-120:-anchor_up]:
													XP_tag = "XP:Z:T"

							new_tag_line = ["AS:i:0", "XO:i:0", "XG:i:0", new_MD_tag, NM_tag, XS_tag, "NH:i:1", XN_tag, XI_tag, XP_tag, XM_tag]
																								
							new_sam = [row[0], new_flag, chr, str(ali_start), row[4], new_cigar] + row[6:9] + [new_seq, new_phred] + new_tag_line										
										
							if ali_flag == "0":
								sense_SJ.write("\t".join(new_sam)+"\n")
#								if len(chr_starts)==0 and MAPQ != "0":
#									sense_SJ_uniq.write("\t".join(new_sam)+"\n") 

							if ali_flag == "16":
								antisense_SJ.write("\t".join(new_sam)+"\n")
#								if len(chr_starts)==0 and MAPQ != "0":
#									antisense_SJ_uniq.write("\t".join(new_sam)+"\n")
								
					except KeyError:
						pass
																											
							
				else:
					chr = tname
								

					XP_tag = "XP:Z:F"	#poliA tag
					XM_tag = "XM:Z:F"	#Adaptador genomico
								
								
					if read+tname+str(start) in genomic_adapter:
						XM_tag = "XM:Z:T"
									
							
					S = start - 1
								
					if ali_flag == "0":
							
						if matches < 95:
							if 'AAAAAA' in Genome[chr][S+matches:S+matches+10]:
								XP_tag = "XP:Z:T"
						else:
							if 'AAAAAA' in Genome[chr][S+matches:S+matches+120]:
								XP_tag = "XP:Z:T"
												
					if ali_flag == "16":									
																												
						if matches < 95:
							if 'TTTTTT' in Genome[chr][S-10:S]:
								XP_tag = "XP:Z:T"
						else:
							if 'TTTTTT' in Genome[chr][S-120:S]:
								XP_tag = "XP:Z:T"
								
								
					new_tag_line = [ XP_tag, XM_tag]
					row = row + new_tag_line			

					genomic.write("\t".join(row)+"\n")

				


#python ~/my_src/Roberto/process_SAM_smallRNAs.py ~/db/genome/hg19.fa wgEncodeCshlShortRnaSeqA549CellShorttotalRawDataRep3.fastq.clip.trim.sam wgEncodeCshlShortRnaSeqA549CellShorttotalRawDataRep3.fastq.clip.trim.sam.adapter Galaxy6-\[Filter_sequences_by_length_on_data_5\].fasta ~/db/transcriptome/hg19/Gene_models/RefSeq/refMrna.fa ~/db/transcriptome/hg19/Gene_models/RefSeq/RefSeq_gen_names ~/db/transcriptome/hg19/Gene_models/RefSeq/RefSeqLinks ~/db/transcriptome/hg19/Gene_models/RefSeq/RefSeq.bed12

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])	
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])

genomic.close()
sense_SJ.close()
antisense_SJ.close()
#sense_SJ_uniq.close()
#antisense_SJ_uniq.close()
