import sys
import csv
import numpy as np
#import matplotlib.pyplot as plt
#from pylab import *

#rcParams['figure.figsize'] = 17, 4.5	

def binP(N, p, x1, x2):
    p = float(p)
    q = p/(1-p)
    k = 0.0
    v = 1.0
    s = 0.0
    tot = 0.0

    while(k<=N):
            tot += v
            if(k >= x1 and k <= x2):
                    s += v
            if(tot > 10**30):
                    s = s/10**30
                    tot = tot/10**30
                    v = v/10**30
            k += 1
            v = v*q*(N+1-k)/k
    return s/tot

def calcBin(vx, vN, vCL = 95):
	'''
	Calculate the exact confidence interval for a binomial proportion

	Usage:
	>>> calcBin(13,100)    
	(0.07107391357421874, 0.21204372406005856)
	>>> calcBin(4,7)   
	(0.18405151367187494, 0.9010086059570312)
	'''


 
	vx = float(vx)
	vN = float(vN)
	#Set the confidence bounds
	vTU = (100 - float(vCL))/2
	vTL = vTU

	vP = vx/vN
	if(vx==0):
		dl = 0.0
	else:
		v = vP/2
		vsL = 0
		vsH = vP
		p = vTL/100

		while((vsH-vsL) > 10**-5):
			if(binP(vN, v, vx, vN) > p):
				vsH = v
				v = (vsL+v)/2
			else:
				vsL = v
				v = (v+vsH)/2
		dl = v

	if(vx==vN):
		ul = 1.0
	else:
		v = (1+vP)/2
		vsL =vP
		vsH = 1
		p = vTU/100
		while((vsH-vsL) > 10**-5):
			if(binP(vN, v, 0, vx) < p):
				vsH = v
				v = (vsL+v)/2
			else:
				vsL = v
				v = (v+vsH)/2
		ul = v
            
	print '%s\t%s\t%s' % (vP, vP - dl, ul - vP)	
#	return vP, vP - dl, ul - vP



def tissiue_bin (means, low, up,tissiue):
	""" Adapta el imput del los ratios, a la funcion para calcular los intervalos de confianza """

	INC_coverage = int(tissiue.split("/")[0])
	alt_coverage = int(tissiue.split("/")[1])	
	intervals = 0, 0, 0

	
	if (INC_coverage + alt_coverage)!=0:
		
		TOTAL = (INC_coverage + alt_coverage)
		intervals = calcBin(INC_coverage,INC_coverage + alt_coverage)
		
	means.append(intervals[0])
	low.append(intervals[1])
	up.append(intervals[2])
	


def main (ratios):
	""" Grafica los ratios de los tejidos """
	

	
	
	reader1 = csv.reader(open(ratios), delimiter = ' ')
	reader1.next()
	
	for row in reader1:
		
		
		tissues_means = []
		tissues_low = []
		tissues_up = []
		
		gene = row[0]
		INC = row[1]
		INC_dn = row[2]
		alt_type = row[3]
		alt_donor = row[4]
		alt_aceptor = row[5]
		alt_intron = row[6]
		alt_dn = row[7]
		
		adipose = row[8]
		adrenal = row[9]
		brain = row[10]
		breast = row[11]
		colon = row[12]
		heart = row[13]
		kidney = row[14]
		liver = row[15]
		lung = row[16]
		lymph_node = row[17]
		ovary = row[18]
		prostate = row[19]
		skeletal_muscle = row[20]
		testes = row[21]
		thyroid = row[22]
		white_blood_cells = row[23]
		
		#if gene=="?":
		#	gene="none"
		
		
		
		if gene=="GKR6" and INC_dn=="GTCC": 
#		if INC_dn!="GTAG" and INC_dn!="GCAG" and INC_dn!="ATAC":	
		
			tissiue_bin(tissues_means, tissues_low, tissues_up, adipose)
			tissiue_bin(tissues_means, tissues_low, tissues_up, adrenal)		
			tissiue_bin(tissues_means, tissues_low, tissues_up, brain)				
			tissiue_bin(tissues_means, tissues_low, tissues_up, breast)	
			tissiue_bin(tissues_means, tissues_low, tissues_up, colon)		
			tissiue_bin(tissues_means, tissues_low, tissues_up, heart)	
			tissiue_bin(tissues_means, tissues_low, tissues_up, kidney)	
			tissiue_bin(tissues_means, tissues_low, tissues_up, liver)	
			tissiue_bin(tissues_means, tissues_low, tissues_up, lung)
			tissiue_bin(tissues_means, tissues_low, tissues_up, lymph_node)
			tissiue_bin(tissues_means, tissues_low, tissues_up, ovary)
			tissiue_bin(tissues_means, tissues_low, tissues_up, prostate)
			tissiue_bin(tissues_means, tissues_low, tissues_up, skeletal_muscle)
			tissiue_bin(tissues_means, tissues_low, tissues_up, testes)
			tissiue_bin(tissues_means, tissues_low, tissues_up, thyroid)
			tissiue_bin(tissues_means, tissues_low, tissues_up, white_blood_cells)		
			

#			max_up = max(tissues_up)
#			
#			ind = np.arange(16)
#			width = 0.7 
#			fig = plt.bar(ind, tissues_means,   width, color='y', yerr=[tissues_low, tissues_up])
#			fig = plt.bar(ind, tissues_means,   width, color='r', yerr=[tissues_low, tissues_up])		

#			Title = " ".join(row[:8])

#			plt.ylabel('(Non-canonical intron)/(Splicing variant)')
#			plt.ylabel('(Non-canonical intron)/(Total splicing variant)')
			
#			plt.title(Title)
#			plt.xticks(ind+width/2, ("Adipose", "Adrenal", "Brain", "Breast", "Colon", "Heart", "Kidney", "Liver", "Lung", "L.Node", "Ovary", "Prostate", "S.Muscle", "Testes", "Thyroid",  "W.B.Cell") )
#			plt.yticks(np.arange(0,1.01,0.1))


			
#			file_name = gene.split(".")[0] + "_" + INC + "_" + alt_intron
			
			names = ["Adipose", "Adrenal", "Brain", "Breast", "Colon", "Heart", "Kidney", "Liver", "Lung", "L.Node", "Ovary", "Prostate", "S.Muscle", "Testes", "Thyroid",  "W.B.Cell"]
			print gene, INC_dn, alt_dn
			for t in [zip(names,tissues_means, tissues_low, tissues_up)[i] for i in (0,3,5,11,12,13,15)]:
				print '\t'.join(map(str, t))			
			



			
			
#			plt.show()
			

			#savefig(file_name)
#			plt.clf()

		
if __name__ == '__main__':
#	main(sys.argv[1])
	calcBin(66,66+108+42)
