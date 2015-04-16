import sys
import csv
from collections import defaultdict
from operator import itemgetter, attrgetter

GENCODE_TSS = set([])


def TSS_loader(gencode_TSS):

	for row in csv.reader(open(gencode_TSS), delimiter = '\t'):

		chrom, start, end, name, score, strand = row

		TSS = chrom + ":" + start 

		GENCODE_TSS.add(TSS)




def main(CIR_TF_Dnase):

	CIR_info = {}
	CIR_names = set([])

	TFs = {}
	Dnase = {}


	for row in csv.reader(open(CIR_TF_Dnase), delimiter = '\t'):

		CIR_chrom = row[0]
		CIR_start = row[1]
		CIR_end = row[2]
		CIR_name = row[3]
		CIR_score = row[4]
		CIR_intersect_type = row[6]
		#CIR_intersect_type_info =  " ".join(row[6:])

		CIR_info[CIR_name] = (CIR_chrom, CIR_start, CIR_end, CIR_score)
		CIR_names.add(CIR_name)

		if CIR_intersect_type == "TFs":

			TFs[CIR_name] = "|".join(row[7:-1])

		if CIR_intersect_type == "Dnase":
			Dnase[CIR_name] = "|".join(row[7:-3])

	Final_Table = []

	

	for CIR_name in CIR_names:

		is_in_promoter = False

		CIR_chrom, CIR_start, CIR_end, CIR_score = CIR_info[CIR_name]

		TFs_info = "."
		Dnase_info = "."

		TF_number = 0

		try:
			TFs_info = TFs[CIR_name]
			TF_number = int(TFs_info.split("|")[-1])

			TF_chrom = TFs_info.split("|")[0]
			TF_start = int(TFs_info.split("|")[1])
			TF_end = int(TFs_info.split("|")[2])

			for n in range(TF_start, TF_end):

				pos = TF_chrom + ":" + str(n)

				if pos in GENCODE_TSS:
					is_in_promoter = True

				# if CIR_name == "18975_ARRDC3":
				# 	if pos
				# 	print TF_chrom, TF_start, TF_end
					#print TF_chrom, TF_start, TF_end, CIR_chrom, CIR_start, CIR_end, CIR_name, CIR_score


		except KeyError:
			pass



		try:
			Dnase_info = Dnase[CIR_name]

		except KeyError:
			pass


		if is_in_promoter==False:

			Final_Table.append([CIR_chrom, CIR_start, CIR_end, CIR_name, CIR_score, TFs_info, Dnase_info, TF_number])

	
	for row in sorted(Final_Table, key=itemgetter(-1) , reverse=True):
		print "\t".join(row[:-1])











if __name__ == '__main__':
	TSS_loader(sys.argv[2])	
	main(sys.argv[1])


# chr12   66289984        66290165        5244_HMGA2      181     .       Dnase   chr12   66289680        66290215        42      915     48      0,9,19,33,82,81,83,92,118,119,4,8,82,16,32,37,36,34,35,38,39,41,42,43,49,61,63,67,67,68,68,69,70,76,96,100,101,103,108,108,111,1,1,31,79,79,71,73,      107,286,413,161,271,365,91,571,254,195,56,99,271,59,450,141,191,619,313,250,98,227,122,134,73,60,92,810,782,698,509,115,655,267,116,171,433,206,259,184,261,56,50,134,915,355,63,72,
# chr12   66289984        66290165        5244_HMGA2      181     .       TFs     chr12   66289298        66290665        HDAC2,YY1,MAFK,RBBP5,ELF1,MYBL2,MAFF,TEAD4,CHD1,CEBPB,ARID3A,REST,POLR2A,EP300,NFIC,MXI1,TBP,TCF7L2,FOXA1,MBD4,FOXA2,TCF12,BCL3 23      .
