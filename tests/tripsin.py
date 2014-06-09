


cuts = [0]
f = []

def main(s):

	for i in range(len(s)-1):
			if ( (s[i]=='K' or s[i]=='R')and (s[i+1]!='P')):
					cuts.append(i)
	cuts.append(len(s))		
					
	for a, b in zip(cuts, cuts[1:]):
		f.append(P[a:b])
	
	print f
			
main("LTRPTGKJHIKPTHHKTTGHV")
