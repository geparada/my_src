import csv
import sys
import numpy as np
from collections import defaultdict
from scipy import stats

import numpy as np
import scipy as sp
from scipy.stats import chi2_contingency
from math import log



def main (contingency_table):
	""" Calcula estadisticas de una tabal de contingencia 2x2 """
	
	SRS_types = set([])
	tables  = {}
	
	for row in csv.reader(open(contingency_table), delimiter = '\t'):
		
		ID, non_can, can = row
		SRS, tag = ID.split("_")
		
		SRS_types.add(SRS)
		tables[ID] = [int(non_can), int(can)]
	
	for srs in SRS_types:
		table = []
		table.append(tables[srs + "_YES"])
		table.append(tables[srs + "_NO"])
		
		
		
		obs = np.array(table)
		chi2, chi2_pvalue, chi2_dof, chi2_ex = chi2_contingency(obs, correction=False)
		chi2_yates, chi2_yates_pvalue, chi2_yates_dof, chi2_yates_ex = chi2_contingency(obs, correction=True)
		fisher_oddsratio, fisher_pvalue = stats.fisher_exact(table)
		
#		print srs, table, fisher_oddsratio, fisher_pvalue, chi2, chi2_pvalue, chi2_dof, chi2_ex
		
		print srs, fisher_oddsratio, log(fisher_oddsratio, 2), fisher_pvalue, chi2, chi2_pvalue, chi2_yates, chi2_yates_pvalue
		 
	


if __name__ == '__main__':
	main(sys.argv[1])	
