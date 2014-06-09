import sys
import csv


def DRcounter(file):
	reader = csv.reader(open(file), dialect='excel-tab' )
	csv.field_size_limit(1000000000)
	for row in reader:
		SJ5=row[15]
		SJ3=row[16]
		L=len(SJ5)/2
		SJ5U = SJ5[:L]
		SJ5D = SJ5[L:]
		SJ3U = SJ3[:L]
		SJ3D = SJ3[L:]
		DRU = 0
		DRD = 0
		if SJ5U==SJ3U:
			DRU = L
		else:
			while SJ5U[L-1-DRU]==SJ3U[L-1-DRU]:
				DRU = (DRU + 1)
				if  SJ5U[L-1-DRU]!=SJ3U[L-1-DRU]: 
					break
		if SJ5D==SJ3D:
			DRD = L
		else:
			 while SJ5D[DRD]==SJ3D[DRD]:
                                DRD = (DRD + 1)
				if SJ5D[DRD]!=SJ3D[DRD]:
					break


		print row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], SJ5U, SJ5D, SJ3U, SJ3D, DRU, DRD, DRU+DRD	
 		#print SJ5U, SJ5D, " "*5 ,SJ3U, SJ3D , DRU, DRD, DRU+DRD


if __name__ == '__main__':                                      #Permite tomar argumentos mediante el terminal
    DRcounter(sys.argv[1])
