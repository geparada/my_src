import sys
import csv

def percent(x, total):
	return 100*float(x)/float(total)

def compare(clip, trim):

	reader1 = csv.reader(open(clip), delimiter = ' ')
	reader2 = csv.reader(open(trim), delimiter = ' ')

	
	introns_clip = []
	clip_list = []	
	for row in reader1:
		read = row[1]
		data = [row[2], row[3]]
		clip_list.append((read, data))
		introns_clip.append(row[1])

	clip = dict(clip_list)

        introns_trim = []
        trim_list = []
        for row in reader2:
                read = row[1]
                data = [row[2], row[3]]
                trim_list.append((read, data))
                introns_trim.append(row[1])

        trim = dict(trim_list)
	
	set_clip = set(introns_clip)
	
	c_clip = 0
	c_clipTP = 0
	c_clipFP = 0
	for c in set_clip:
		c_clip += 1
		if clip[c][0] == 'Known':
			c_clipTP += 1			
		elif clip[c][1] != 'GTAG' and clip[c][1] != 'GCAG' and clip[c][1] != 'ATAC':
			c_clipFP += 1

	print 'Dataset', 'Total', 'Canonical&Known', 'Non-Canonical&Unknown', 'Canonical&Unknown' 
	print 'Clip', c_clip, percent(c_clipTP, c_clip) , percent (c_clipFP, c_clip), percent(c_clip - c_clipTP - c_clipFP, c_clip)


        set_trim = set(introns_trim)

        c_trim = 0
        c_trimTP = 0
        c_trimFP = 0
        for t in set_trim:
                c_trim += 1

                if trim[t][0] == 'Known':
                        c_trimTP += 1
                elif trim[t][1] != 'GTAG' and trim[t][1] != 'GCAG' and trim[t][1] != 'ATAC':
                        c_trimFP += 1

        print 'Trim', c_trim, percent(c_trimTP, c_trim), percent(c_trimFP, c_trim),  percent(c_trim - c_trimTP - c_trimFP, c_trim)


	only_clip = set(introns_clip) - set(introns_trim)
	
	c_clip = 0
	c_clipTP = 0
	c_clipFP = 0
	for c in only_clip:
		c_clip += 1
		if clip[c][0] == 'Known':
			c_clipTP += 1			
		elif clip[c][1] != 'GTAG' and clip[c][1] != 'GCAG' and clip[c][1] != 'ATAC':
			c_clipFP += 1

	print 'Solo_Clip', c_clip, percent(c_clipTP, c_clip) , percent (c_clipFP, c_clip), percent(c_clip - c_clipTP - c_clipFP, c_clip)


        only_trim = set(introns_trim) - set(introns_clip)

        c_trim = 0
        c_trimTP = 0
        c_trimFP = 0
        for t in only_trim:
                c_trim += 1

                if trim[t][0] == 'Known':
                        c_trimTP += 1
                elif trim[t][1] != 'GTAG' and trim[t][1] != 'GCAG' and trim[t][1] != 'ATAC':
                        c_trimFP += 1

        print 'Solo_Trim', c_trim, percent(c_trimTP, c_trim), percent(c_trimFP, c_trim),  percent(c_trim - c_trimTP - c_trimFP, c_trim)






if __name__ == '__main__':
	compare(sys.argv[1], sys.argv[2])
