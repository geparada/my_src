import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def main (mast_line, chr, start, end, l):
	start = int(start)
	l = int(l)
	n = 0
	c = 0
	for m in mast_line.split("]_"):

		try:
			i, strand = m.split("_[")
			i = int(i)
			strand = strand[0]
			n += i
			correction_factor = l*c 
			c += 1			
			name = "Motif_" + str(c)			
			print '%s\t%s\t%s\t%s\t%s\t%s' % (chr, start + n + correction_factor, start + n + l + correction_factor, name, "0", strand)
			
		except ValueError:
			pass 



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
	#main('29044_[+1]_4925_[+1]_1504_[+1]_3259_[-1]_5641_[-1]_7387_[+1]_5760_[-1]_9970_[+1]_1701_[+1]_597_[-1]_13109_[+1]_2969_[-1]_7205_[+1]_6708', 'chr8', '40666824', '40766824', '17')




#python ~/my_src/Fran/mast_to_genome.py 29044_[+1]_4925_[+1]_1504_[+1]_3259_[-1]_5641_[-1]_7387_[+1]_5760_[-1]_9970_[+1]_1701_[+1]_597_[-1]_13109_[+1]_2969_[-1]_7205_[+1]_6708 chr8 40666824 40766824 17

#extra 99_[+1]_10237_[-1]_175_[-1]_2627_[-1]_6914_[-1]_12863 chr8 40760000 40793000 17


#29044_[+1]_4925_[+1]_1504_[+1]_3259_[-1]_5641_[-1]_7387_[+1]_5760_[-1]_9970_[+1]_1701_[+1]_597_[-1]_13109_[+1]_2969_[-1]_7205_[+1]_6708
