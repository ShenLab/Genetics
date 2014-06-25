import csv

#read  data
f = open('PCGC_summary_readable.csv','Ur')
r = csv.reader(f)
f1 = open('PCGC_summary.csv','wb')
w = csv.writer(f1)
head = r.next()
w.writerow(head)

for line in r:
    if int(line[0]) <= 35:# and line[16] in ['D','NA'] and line[18] in ['D','P','NA']:
        w.writerow(line)
    
f1.close()
