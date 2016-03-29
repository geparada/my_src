import sys
import csv
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collection import defaultdict

#dicionario de fag [Rd1/Rd2, +/-]
flag_dict = {'73':[1,1], '89':[1,1], '121':[1,-1], '153':[-1,-1], '185':[-1,-1], '137':[-1,1], '99':[1,1], '147':[-1,-1], '83':[1,-1], '163':[-1,1], '67':[1,1], '115':[1,-1], '179':[-1,-1], '81':[1,-1], "161":[-1,1], '97':[1,1], '145':[-1,-1], '65':[1,1], '129':[-1,1], '113':[1,-1], '177':[-1,-1] }

Genome = {}
non_can = set([])
can = set([])

dn_count = defaultdict(int)

def Genomictabulator(fasta):
	
	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",	

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq
		
	print >> sys.stderr, "OK"

	f.close()


def main(bam, forward):

	samfile = pysam.Samfile(bam, "rb")
	
	#print samfile.gettid("0")

	for read in samfile.fetch():
		start = read.pos
		flag = read.flag
		
		
		chr = samfile.getrname(read.tid)		

		Exon_starts = [start]
		Exon_ends = []
				
		block = 0
		var_index = 0
		
		pair_ori = 0
		if forward == "Rd1":
			pair_ori = 1
		elif forward == "Rd2":
			pair_ori = -1
						
		self_strand = 1
		pair_strand = '+'
							
		#Si no se tiene el flag XS:A:- se tienen que implementar las operaciones a nivel de bits. Esto es necesario porque los intrones no canonicos no tienen ese tag.
				
		if (1 & int(flag)):    #paired end
			pair_number = flag_dict[str(flag)][0]
			self_strand = flag_dict[str(flag)][1]
			if pair_ori*self_strand*pair_number==-1:
				pair_strand = '-'
											
		elif (16 & int(flag)):   #single end      
			self_strand = -1
			pair_strand = '-'
								
	
		for cigar_tuple in read.cigar:
			
				
			var_type = cigar_tuple[0]    #traduciendo de Bam a Op
			if cigar_tuple[0] == 0:
				var_type = 'M'
			elif cigar_tuple[0] == 1:
				var_type = 'I'				
			elif cigar_tuple[0] == 2:
				var_type = 'D'				
			elif cigar_tuple[0] == 3:
				var_type = 'N'	
			
			var_value = cigar_tuple[1]
			var_index += 1
					
			if var_type == 'M':
				block += var_value						
						
			if var_type == 'D':
				block += var_value
											
			if var_type == 'I':
				block += 0
						
			if var_type == 'N':
				Exon_ends.append(Exon_starts[-1] + block)
				Exon_starts.append(Exon_ends[-1] + var_value)
				block = 0
						
			if var_index == len(cigar_tuple):
				Exon_ends.append(Exon_starts[-1] + block)

			for e5s, e5e, e3s, e3e  in zip(Exon_starts, Exon_ends, Exon_starts[1:], Exon_ends[1:]):
				e5len= e5e - e5s 
				e3len = e3e - e3s
				istart = e5e
				iend = e3s
				ilen = iend - istart
				intron = chr + ":" +  str(istart) + pair_strand + str(iend)

				dn = Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]

				if pair_strand == '-':
					dn = dn.reverse_complement()
							
				dn = str(dn).upper()
				
				if dn=="GTAG" or dn=="GCAG" or dn=="ATAC":
					can.add(intron)
				else:
					non_can.add(intron)

				dn_count[dn]+=1
	
	print "canonical introns", len(can)
	print "non-canonical introns", len(non_can)

	dn_table = dn_count.items()
	dn_table.sort(key=lambda x: x[1])

	for i in reversed(dn_table):
		dn = i[0]
		N = i[1]
		print dn, N

	samfile.close()



if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3])
