from collections import defaultdict
import sys

def uniq(file):
	data = defaultdict(list)
	with open(file, 'rb') as f:
	    for line in sorted(f.readlines()):
		data[line[0]].append(line)
	for key in sorted(data.iterkeys()):
	    if len(data[key]) == 1:
		print data[key]

if __name__ == "__main__":
    uniq(sys.argv[1])
